#include <algorithm>

#include "SDL2/SDL.h"

#include "../../lib/triangle/triangle.h"
#include "../triangle/util.h"

#include <constants.hpp>

#include "mesh.hpp"

namespace nikfemm {
    Mesh::Mesh() {
        
    }

    Mesh::~Mesh() {

    }

    void Mesh::plot() {
        // draw the mesh
        // returns zero on success else non-zero
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            printf("error initializing SDL: %s\n", SDL_GetError());
        }
        SDL_Window* win = SDL_CreateWindow("GAME", // creates a window
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        2000, 2000, 0);
            
        // triggers the program that controls
        // your graphics hardware and sets flags
        Uint32 render_flags = SDL_RENDERER_ACCELERATED;
    
        // creates a renderer to render our images
        SDL_Renderer* rend = SDL_CreateRenderer(win, -1, 0);

        // clears the window
        SDL_RenderClear(rend);

        // get mesh enclosing rectangle
        double min_x = 1000000;
        double min_y = 1000000;
        double max_x = -1000000;
        double max_y = -1000000;
        for (auto v : vertices) {
            if (v->p.x < min_x) {
                min_x = v->p.x;
            }
            if (v->p.y < min_y) {
                min_y = v->p.y;
            }
            if (v->p.x > max_x) {
                max_x = v->p.x;
            }
            if (v->p.y > max_y) {
                max_y = v->p.y;
            }
        }

        // object to window ratio
        double ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        double x_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        double y_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        double x_offset = 0.5 * 2000 - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        double y_offset = 0.5 * 2000 - 0.5 * (max_y + min_y) * y_scale;

        // render

        while(true){
            // draw the triangles
            for (auto e : elements) {
                SDL_SetRenderDrawColor(rend, 255, 0, 0, 255);

                double p1x = e->vertices[0]->p.x;
                double p1y = e->vertices[0]->p.y;
                double p2x = e->vertices[1]->p.x;
                double p2y = e->vertices[1]->p.y;
                double p3x = e->vertices[2]->p.x;
                double p3y = e->vertices[2]->p.y;

                SDL_RenderDrawLine(rend, x_offset + p1x * x_scale, y_offset + p1y * y_scale, x_offset + p2x * x_scale, y_offset + p2y * y_scale);
                SDL_RenderDrawLine(rend, x_offset + p2x * x_scale, y_offset + p2y * y_scale, x_offset + p3x * x_scale, y_offset + p3y * y_scale);
                SDL_RenderDrawLine(rend, x_offset + p3x * x_scale, y_offset + p3y * y_scale, x_offset + p1x * x_scale, y_offset + p1y * y_scale);
            }

            // draw the points
            for (auto v : vertices) {
                SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
                // draw a square centered at the point
                SDL_Rect rect;
                rect.x = x_offset + v->p.x * x_scale - 2;
                rect.y = y_offset + v->p.y * y_scale - 2;
                rect.w = 4;
                rect.h = 4;
                SDL_RenderFillRect(rend, &rect);
            }

            SDL_RenderPresent(rend);

            SDL_Event event;
            if (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    break;
                }
            }

            // sleep for 1 second
            SDL_Delay(10);
        }
    }

    void Mesh::mesh(Drawing &drawing) {
        triangulateio in, out;

        /* set up the input */
        struct TempSegment {
            uint64_t id1;
            uint64_t id2;
        };

        in.numberofpoints = drawing.segments.size() * 2;
        in.pointlist = (TRI_REAL*)malloc(in.numberofpoints * 2 * sizeof(TRI_REAL));
        in.numberofsegments = drawing.segments.size();
        in.segmentlist = (int*)malloc(in.numberofsegments * 2 * sizeof(int));

        uint32_t points_idx = 0;
        uint32_t segments_idx = 0;

        // for each segment, check if one of the endpoints is already in the list, if not, add it
        for (auto s : drawing.segments) {
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
        in.pointmarkerlist = NULL;
        in.numberofpointattributes = 0;
        in.pointattributelist = NULL;
        in.segmentmarkerlist = NULL;
        in.numberofholes = 0;
        in.holelist = NULL;

        in.numberofregions = drawing.regions.size();
        in.regionlist = (TRI_REAL*)malloc(in.numberofregions * 4 * sizeof(TRI_REAL));
        uint64_t i = 0;
        for (auto r : drawing.regions) {
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

        triangulate(switches, &in, &out, NULL);

        /* copy the triangles */
        // for each vertex in each triangle, add a point to the list
        for (uint64_t i = 0; i < out.numberoftriangles; i++) {
            uint64_t p1 = out.trianglelist[3 * i];
            uint64_t p2 = out.trianglelist[3 * i + 1];
            uint64_t p3 = out.trianglelist[3 * i + 2];

            TriangleVertex* v1 = new TriangleVertex(out.pointlist[2 * p1], out.pointlist[2 * p1 + 1]);
            TriangleVertex* v2 = new TriangleVertex(out.pointlist[2 * p2], out.pointlist[2 * p2 + 1]);
            TriangleVertex* v3 = new TriangleVertex(out.pointlist[2 * p3], out.pointlist[2 * p3 + 1]);

            std::pair<std::unordered_set<TriangleVertex*>::iterator, bool> ret = vertices.insert(v1);
            boundary_vertices.reserve(BOUNDARY_VERTICES);
            if (out.pointmarkerlist[p1] == 1) {
                boundary_vertices.push_back(v1);
            }
            if (!ret.second) {
                delete v1;
                v1 = *ret.first;
            }
            ret = vertices.insert(v2);
            if (out.pointmarkerlist[p2] == 1) {
                boundary_vertices.push_back(v2);
            }
            if (!ret.second) {
                delete v2;
                v2 = *ret.first;
            }
            ret = vertices.insert(v3);
            if (out.pointmarkerlist[p3] == 1) {
                boundary_vertices.push_back(v3);
            }
            if (!ret.second) {
                delete v3;
                v3 = *ret.first;
            }

            // add the triangle
            TriangleElement* t = new TriangleElement(v1, v2, v3);
            elements.insert(t);

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

        // sort boundary vertices
        auto comparator = TriangleVertex::atanCompare(center);
        std::sort(boundary_vertices.begin(), boundary_vertices.end(), comparator);
    }

    void Mesh::addKelvinBoundaryConditions() {
        
    }
}