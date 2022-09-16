#include <unordered_set>
#include <set>
#include <chrono>
#include <array>

#include <constants.hpp>

#include "../lib/triangle/triangle.h"
#include "triangle/util.h"

#include "drawing.hpp"
#include "mesh/mesh.hpp"

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

        in.numberofpoints = segments.size() * 2;
        in.pointlist = (TRI_REAL*)malloc(in.numberofpoints * 2 * sizeof(TRI_REAL));
        in.numberofsegments = segments.size();
        in.segmentlist = (int*)malloc(in.numberofsegments * 2 * sizeof(int));

        uint32_t points_idx = 0;
        uint32_t segments_idx = 0;

        // for each segment, check if one of the endpoints is already in the list, if not, add it
        for (auto s : segments) {
            bool p1_found = false;
            bool p2_found = false;
            uint32_t p1_id = 0;
            uint32_t p2_id = 0;
            uint32_t i = 0;
            for (uint32_t i = 0; i < points_idx; i++) {
                // compare using epsilon
                if (fabs(in.pointlist[i * 2] - s.p1.x) < EPSILON && fabs(in.pointlist[i * 2 + 1] - s.p1.y) < EPSILON) {
                    p1_found = true;
                    p1_id = i;
                }
                if (fabs(in.pointlist[i * 2] - s.p2.x) < EPSILON && fabs(in.pointlist[i * 2 + 1] - s.p2.y) < EPSILON) {
                    p2_found = true;
                    p2_id = i;
                }
                if (p1_found && p2_found) {
                    break;
                }
            }
            if (!p1_found) {
                in.pointlist[points_idx * 2] = s.p1.x;
                in.pointlist[points_idx * 2 + 1] = s.p1.y;
                p1_id = points_idx;
                points_idx++;
            }
            if (!p2_found) {
                in.pointlist[points_idx * 2] = s.p2.x;
                in.pointlist[points_idx * 2 + 1] = s.p2.y;
                p2_id = points_idx;
                points_idx++;
            }
            // add the segment
            in.segmentlist[segments_idx * 2] = p1_id;
            in.segmentlist[segments_idx * 2 + 1] = p2_id;
            segments_idx++;
        }

        in.numberofpoints = points_idx + 1;

        auto start1 = std::chrono::high_resolution_clock::now();

        in.pointmarkerlist = NULL;
        in.numberofpointattributes = 0;
        in.pointattributelist = NULL;
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

        char switches[] = "pzq20cAQa0.01";

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
        for (uint64_t i = 0; i < out.numberoftriangles; i++) {
            uint64_t p1 = out.trianglelist[3 * i];
            uint64_t p2 = out.trianglelist[3 * i + 1];
            uint64_t p3 = out.trianglelist[3 * i + 2];

            TriangleVertex* v1 = new TriangleVertex(out.pointlist[2 * p1], out.pointlist[2 * p1 + 1]);
            TriangleVertex* v2 = new TriangleVertex(out.pointlist[2 * p2], out.pointlist[2 * p2 + 1]);
            TriangleVertex* v3 = new TriangleVertex(out.pointlist[2 * p3], out.pointlist[2 * p3 + 1]);
            auto astart3 = std::chrono::high_resolution_clock::now();

            std::pair<std::unordered_set<TriangleVertex*>::iterator, bool> ret = mesh->vertices.insert(v1);
            if (out.pointmarkerlist[p1] == 1) {
                mesh->boundary_vertices.insert(v1);
            }
            if (!ret.second) {
                delete v1;
                v1 = *ret.first;
            }
            ret = mesh->vertices.insert(v2);
            if (out.pointmarkerlist[p2] == 1) {
                mesh->boundary_vertices.insert(v2);
            }
            if (!ret.second) {
                delete v2;
                v2 = *ret.first;
            }
            ret = mesh->vertices.insert(v3);
            if (out.pointmarkerlist[p3] == 1) {
                mesh->boundary_vertices.insert(v3);
            }
            if (!ret.second) {
                delete v3;
                v3 = *ret.first;
            }

            // add the triangle
            TriangleElement* t = new TriangleElement(v1, v2, v3);
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

        return mesh;
    }
}