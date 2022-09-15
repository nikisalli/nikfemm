#include <unordered_set>
#include <set>
#include <chrono>

#include <constants.hpp>

#include "../lib/triangle/triangle.h"
#include "triangle/util.h"

#include "drawing.hpp"
#include "mesh.hpp"

namespace nikfemm {
    Drawing::Drawing() {
        
    }

    Drawing::~Drawing() {
        
    }

    void Drawing::drawRectangle(Point p1, Point p2) {
        segments.insert(Segment(p1, Point(p2.x, p1.y)));
        segments.insert(Segment(Point(p2.x, p1.y), p2));
        segments.insert(Segment(p2, Point(p1.x, p2.y)));
        segments.insert(Segment(Point(p1.x, p2.y), p1));
    }

    void Drawing::drawRectangle(Point p1, double width, double height) {
        drawRectangle(p1, Point(p1.x + width, p1.y + height));
    }

    void Drawing::drawCircle(Point p, double radius, uint32_t n_segments) {
        double angle = 0;
        double angle_step = 2 * PI / n_segments;
        Point p1 = Point(p.x + radius, p.y);
        Point p2;
        for (uint32_t i = 0; i < n_segments; i++) {
            angle += angle_step;
            p2 = Point(p.x + radius * cos(angle), p.y + radius * sin(angle));
            segments.insert(Segment(p1, p2));
            p1 = p2;
        }
    }

    void Drawing::drawCircle(Circle c, uint32_t n_segments) {
        drawCircle(c.center, c.radius, n_segments);
    }

    void Drawing::drawPolygon(Point* points, uint32_t n_points) {
        // check if the polygon self-intersects
        for (uint32_t i = 0; i < n_points; i++) {
            for (uint32_t j = i + 2; j < n_points; j++) {
                if (i == 0 && j == n_points - 1) {
                    continue;
                }
                if (Segment::segmentsIntersect(points[i], points[(i + 1) % n_points], points[j], points[(j + 1) % n_points])) {
                    printf("Error: polygon self-intersects\n");
                    throw std::invalid_argument("polygon self-intersects");
                }
            }
        }

        for (uint32_t i = 0; i < n_points - 1; i++) {
            segments.insert(Segment(points[i], points[i + 1]));
        }
        segments.insert(Segment(points[n_points - 1], points[0]));
    }

    void Drawing::drawPolyLine(Point* points, uint32_t n_points) {
        for (uint32_t i = 0; i < n_points - 1; i++) {
            segments.insert(Segment(points[i], points[i + 1]));
        }
    }

    void Drawing::drawRegion(Point p, uint32_t region_id) {
        regions.insert(DrawingRegion(p, region_id));
    }

    void Drawing::drawRegion(Point p, PredefinedRegion region) {
        regions.insert(DrawingRegion(p, region.region_id));
    }

    Mesh* Drawing::mesh() {
        auto start = std::chrono::high_resolution_clock::now();

        triangulateio in, out;

        /* set up the input */
        struct TempSegment {
            uint64_t id1;
            uint64_t id2;
        };

        std::vector<Point> points;
        std::vector<TempSegment> temp_segments;
        // for each segment, check if one of the endpoints is already in the list, if not, add it
        printf("segments size: %lu\n", segments.size());
        for (auto s : segments) {
            bool p1_found = false;
            bool p2_found = false;
            uint64_t p1_id = 0;
            uint64_t p2_id = 0;
            uint64_t i = 0;
            for (uint32_t i = 0; i < points.size(); i++) {
                if (points[i] == s.p1) {
                    p1_found = true;
                    p1_id = i;
                }
                if (points[i] == s.p2) {
                    p2_found = true;
                    p2_id = i;
                }
                if (p1_found && p2_found) {
                    break;
                }
            }
            if (!p1_found) {
                points.push_back(s.p1);
                p1_id = points.size() - 1;
            }
            if (!p2_found) {
                points.push_back(s.p2);
                p2_id = points.size() - 1;
            }
            temp_segments.push_back(TempSegment{ p1_id, p2_id });
        }

        auto start1 = std::chrono::high_resolution_clock::now();

        in.numberofpoints = points.size();
        in.pointlist = (TRI_REAL*)malloc(in.numberofpoints * 2 * sizeof(TRI_REAL));
        for (uint64_t i = 0; i < in.numberofpoints; i++) {
            in.pointlist[2 * i] = points[i].x;
            in.pointlist[2 * i + 1] = points[i].y;
        }
        in.pointmarkerlist = NULL;
        in.numberofpointattributes = 0;
        in.pointattributelist = NULL;

        in.numberofsegments = segments.size();
        in.segmentlist = (int*)malloc(in.numberofsegments * 2 * sizeof(int));
        for (uint64_t i = 0; i < in.numberofsegments; i++) {
            in.segmentlist[2 * i] = temp_segments[i].id1;
            in.segmentlist[2 * i + 1] = temp_segments[i].id2;
        }
        in.segmentmarkerlist = NULL;

        in.numberofholes = 0;
        in.holelist = NULL;

        in.numberofregions = regions.size();
        in.regionlist = (TRI_REAL*)malloc(in.numberofregions * 4 * sizeof(TRI_REAL));
        uint64_t i = 0;
        for (auto r : regions) {
            in.regionlist[4 * i] = r.p.x;
            in.regionlist[4 * i + 1] = r.p.y;
            in.regionlist[4 * i + 2] = r.region_id;
            in.regionlist[4 * i + 3] = 0;
            i++;
        }

        out.pointlist = (TRI_REAL*)NULL;

        char switches[] = "pzqceAa0.01";

        /* Make necessary initializations so that Triangle can return a */
        /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

        out.pointlist = (TRI_REAL *) NULL;            /* Not needed if -N switch used. */
        /* Not needed if -N switch used or number of point attributes is zero: */
        out.pointattributelist = (TRI_REAL *) NULL;
        out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
        out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
        /* Not needed if -E switch used or number of triangle attributes is zero: */
        out.triangleattributelist = (TRI_REAL *) NULL;
        out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
        /* Needed only if segments are output (-p or -c) and -P not used: */
        out.segmentlist = (int *) NULL;
        /* Needed only if segments are output (-p or -c) and -P and -B not used: */
        out.segmentmarkerlist = (int *) NULL;
        out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
        out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

        auto start2 = std::chrono::high_resolution_clock::now();

        triangulate(switches, &in, &out, NULL);

        auto start3 = std::chrono::high_resolution_clock::now();

        Mesh* mesh = new Mesh();

        /* copy the triangles */
        // for each vertex in each triangle, add a point to the list
        auto total_time1 = std::chrono::duration_cast<std::chrono::microseconds>(start1 - start).count();
        for (uint64_t i = 0; i < out.numberoftriangles; i++) {
            uint64_t p1 = out.trianglelist[3 * i];
            uint64_t p2 = out.trianglelist[3 * i + 1];
            uint64_t p3 = out.trianglelist[3 * i + 2];

            Vertex* v1 = new Vertex(out.pointlist[2 * p1], out.pointlist[2 * p1 + 1]);
            Vertex* v2 = new Vertex(out.pointlist[2 * p2], out.pointlist[2 * p2 + 1]);
            Vertex* v3 = new Vertex(out.pointlist[2 * p3], out.pointlist[2 * p3 + 1]);

            std::pair<std::unordered_set<Vertex*>::iterator, bool> ret = mesh->vertices.insert(v1);
            if (!ret.second) {
                delete v1;
                v1 = *ret.first;
            }
            ret = mesh->vertices.insert(v2);
            if (!ret.second) {
                delete v2;
                v2 = *ret.first;
            }
            ret = mesh->vertices.insert(v3);
            if (!ret.second) {
                delete v3;
                v3 = *ret.first;
            }

            // add the triangle
            Element* t = new Element(v1, v2, v3);
            mesh->elements.insert(t);

            // add the triangle as a neighbor of each vertex
            v1->addAdjacentElement(t);
            v2->addAdjacentElement(t);
            v3->addAdjacentElement(t);

            // add vertices as neighbors of each other if they are not already
            v1->addAdjacentVertex(v2);
            v1->addAdjacentVertex(v3);
            v2->addAdjacentVertex(v1);
            v2->addAdjacentVertex(v3);
            v3->addAdjacentVertex(v1);
            v3->addAdjacentVertex(v2);
        }
        
        auto start4 = std::chrono::high_resolution_clock::now();

        // add neighbor elements to each element
        for (auto elem : mesh->elements) {
            // for each pair of vertices in the element find the element that shares those vertices
            uint64_t vert_num = elem->vertices.size();
            for (uint64_t j = 0; j < vert_num; j++) {
                for (uint64_t k = j + 1; k < vert_num; k++) {
                    Vertex* v1 = elem->vertices[j];
                    Vertex* v2 = elem->vertices[k];
                    for (uint64_t l = 0; l < v1->adjele.size(); l++) {
                        for (uint64_t m = 0; m < v2->adjele.size(); m++) {
                            if (v1->adjele[l] == v2->adjele[m] && v1->adjele[l] != elem) {
                                elem->addAdjacentElement(v1->adjele[l]);
                            }
                        }
                    }
                }
            }
        }

        auto start5 = std::chrono::high_resolution_clock::now();

        // find the boundary elements
        for (auto elem : mesh->elements) {
            if (elem->adjele.size() < 3) {
                mesh->boundary_elements.insert(elem);
            }
        }

        auto start6 = std::chrono::high_resolution_clock::now();

        // find the boundary vertices
        for (auto elem : mesh->boundary_elements) {
            // for each edge in the element find the edge that is not shared by any other element
            uint64_t vert_num = elem->vertices.size();
            for (uint64_t j = 0; j < vert_num; j++) {
                for (uint64_t k = j + 1; k < vert_num; k++) {
                    Vertex* v1 = elem->vertices[j];
                    Vertex* v2 = elem->vertices[k];
                    bool found = false;
                    for (uint64_t l = 0; l < v1->adjele.size(); l++) {
                        for (uint64_t m = 0; m < v2->adjele.size(); m++) {
                            if (v1->adjele[l] == v2->adjele[m] && v1->adjele[l] != elem) {
                                found = true;
                            }
                        }
                    }
                    if (!found) {
                        mesh->boundary_vertices.insert(v1);
                        mesh->boundary_vertices.insert(v2);
                    }
                }
            }
        }

        auto start7 = std::chrono::high_resolution_clock::now();

        printf("Create temp segments: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start1 - start).count()*1000);
        printf("Fill triangulateio: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start2 - start1).count()*1000);
        printf("Triangulate: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start3 - start2).count()*1000);
        printf("Fill mesh: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start4 - start3).count()*1000);
        printf("Find neighbors e-e: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start5 - start4).count()*1000);
        printf("Find boundary elements: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start6 - start5).count()*1000);
        printf("Find boundary vertices: %f\n", std::chrono::duration_cast<std::chrono::duration<double>>(start7 - start6).count()*1000);

        return mesh;
    }
}