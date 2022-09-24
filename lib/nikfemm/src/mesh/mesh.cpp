#include <algorithm>
#include <stdexcept>
#include <chrono>

#include "SDL2/SDL.h"
#include "opencv2/opencv.hpp"

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
            // clears the window
            SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
            SDL_RenderClear(rend);
            // draw point at x: -1.389164 y: 0.264997
            SDL_SetRenderDrawColor(rend, 0, 255, 0, 255);
            SDL_Rect rect;
            rect.x = x_offset + x_scale * -1.389164;
            rect.y = y_offset + y_scale * 0.264997;
            rect.w = 5;
            rect.h = 5;
            SDL_RenderFillRect(rend, &rect);
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

                for (uint16_t i = 0; i < v->adjvert_count; i++) {
                    SDL_SetRenderDrawColor(rend, 255, 0, 0, 255);
                    SDL_RenderDrawLine(rend, x_offset + v->p.x * x_scale, y_offset + v->p.y * y_scale, x_offset + v->adjvert[i]->p.x * x_scale, y_offset + v->adjvert[i]->p.y * y_scale);
                }
            }

            SDL_RenderPresent(rend);

            SDL_Event event;
            if (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    // destroy the window
                    SDL_DestroyWindow(win);
                    // clean up
                    SDL_Quit();
                    break;
                }
            }

            SDL_KeyboardEvent key = event.key;
            if (key.keysym.sym == SDLK_UP) {
                y_offset += 10;
            } else if (key.keysym.sym == SDLK_DOWN) {
                y_offset -= 10;
            } else if (key.keysym.sym == SDLK_LEFT) {
                x_offset += 10;
            } else if (key.keysym.sym == SDLK_RIGHT) {
                x_offset -= 10;
            } else if (key.keysym.sym == SDLK_w) {
                y_offset += 10;
            } else if (key.keysym.sym == SDLK_s) {
                y_offset -= 10;
            } else if (key.keysym.sym == SDLK_a) {
                y_offset += 10;
            } else if (key.keysym.sym == SDLK_d) {
                y_offset -= 10;
            } else if (key.keysym.sym == SDLK_q) {
                x_scale *= 0.9;
                y_scale *= 0.9;
            } else if (key.keysym.sym == SDLK_e) {
                x_scale *= 1.1;
                y_scale *= 1.1;
            }

            // sleep for 1 second
        }
    }

    void Mesh::mesh(Drawing &drawing) {
        triangulateio in, out;

        /* set up the input */
        struct TempSegment {
            uint64_t id1;
            uint64_t id2;
        };

        in.numberofpoints = drawing.points.size();
        in.pointlist = (TRI_REAL*)malloc(in.numberofpoints * 2 * sizeof(TRI_REAL));
        in.numberofsegments = drawing.segments.size();
        in.segmentlist = (int*)malloc(in.numberofsegments * 2 * sizeof(int));

        for (uint64_t i = 0; i < drawing.points.size(); i++) {
            in.pointlist[2 * i] = drawing.points[i].x;
            in.pointlist[2 * i + 1] = drawing.points[i].y;
        }

        uint64_t i = 0;
        for (auto s : drawing.segments) {
            in.segmentlist[2 * i] = s.p1;
            in.segmentlist[2 * i + 1] = s.p2;
            i++;
        }

        in.pointmarkerlist = NULL;
        in.numberofpointattributes = 0;
        in.pointattributelist = NULL;
        in.segmentmarkerlist = NULL;
        in.numberofholes = 0;
        in.holelist = NULL;

        in.numberofregions = drawing.regions.size();
        in.regionlist = (TRI_REAL*)malloc(in.numberofregions * 4 * sizeof(TRI_REAL));
        i = 0;
        for (auto r : drawing.regions) {
            in.regionlist[4 * i] = r.p.x;
            in.regionlist[4 * i + 1] = r.p.y;
            in.regionlist[4 * i + 2] = r.region_attribute;
            in.regionlist[4 * i + 3] = 0;
            i++;
        }

        out.pointlist = (TRI_REAL*)NULL;

        char switches[] = "pzq20AQa0.01";

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

        vertices.reserve(out.numberofpoints);
        boundary_vertices.reserve(BOUNDARY_VERTICES);
        for (uint64_t i = 0; i < out.numberofpoints; i++) {
            Vertex* v = new Vertex(out.pointlist[2 * i], out.pointlist[2 * i + 1]);
            vertices.push_back(v);
            if (out.pointmarkerlist[i] == 1) {
                boundary_vertices.push_back(v);
            }
        }
        // for each vertex in each triangle, add a point to the list
        for (uint64_t i = 0; i < out.numberoftriangles; i++) {
            uint64_t p1 = out.trianglelist[3 * i];
            uint64_t p2 = out.trianglelist[3 * i + 1];
            uint64_t p3 = out.trianglelist[3 * i + 2];

            Vertex* v1 = vertices[p1];
            Vertex* v2 = vertices[p2];
            Vertex* v3 = vertices[p3];

            // add vertices as neighbors of each other if they are not already
            v1->addAdjacentVertex(v2);
            v1->addAdjacentVertex(v3);
            v2->addAdjacentVertex(v1);
            v2->addAdjacentVertex(v3);
            v3->addAdjacentVertex(v1);
            v3->addAdjacentVertex(v2);

            // add mu_r
            v1->addAdjacentMu(out.triangleattributelist[i]);
            v2->addAdjacentMu(out.triangleattributelist[i]);
            v3->addAdjacentMu(out.triangleattributelist[i]);
        }

        // sort boundary vertices
        auto comparator = Vertex::atanCompare(center);
        std::sort(boundary_vertices.begin(), boundary_vertices.end(), comparator);
        printf("Number of vertices: %lu, boundary vertices: %lu\n", vertices.size(), boundary_vertices.size());
    }

    void Mesh::kelvinTransformCentered() {
        // kelvin transform each vertex in the mesh
        /*

        x* = (R^2 / |x|^2) * x

        */

        double R_squared = radius * radius;
        // printf("R^2 = %f\n", R_squared);

        for (auto it = vertices.begin(); it != vertices.end(); it++) {
            Vertex* v = *it;
            double mag_squared = v->p.x * v->p.x + v->p.y * v->p.y;
            double scale = R_squared / mag_squared;
            // printf("v = (%f, %f) -> (%f, %f), mag = %f, scale = %f\n", v->p.x, v->p.y, v->p.x * scale, v->p.y * scale, mag_squared, scale);
            v->p = v->p * scale;
        }
    }

    void Mesh::addKelvinBoundaryConditions() {
        // kelvin mesh = km
        // this mesh = tm
        // kelvin mesh boundary vertices = kmb
        // this mesh boundary vertices = tmb
        // kelvin mesh vertices = kmv
        // this mesh vertices = tmv

        Mesh kelvin_mesh;
        Drawing kelvin_drawing;

        // for each consecutive pair of boundary vertices, add a segment
        kelvin_drawing.drawSegment(*(boundary_vertices.end() - 1), *(boundary_vertices.begin()));
        for (uint64_t i = 0; i < boundary_vertices.size() - 1; i++) {
            kelvin_drawing.drawSegment(*(boundary_vertices.begin() + i), *(boundary_vertices.begin() + i + 1));
        }
        kelvin_drawing.drawRegion(Point(0, 0), BOUNDARY_REGION);

        kelvin_mesh.mesh(kelvin_drawing);
        kelvin_mesh.center = Point(0, 0);
        kelvin_mesh.radius = radius;

        // merge km into tm
        if (boundary_vertices.size() != kelvin_mesh.boundary_vertices.size()) {
            throw std::runtime_error("Kelvin mesh boundary vertices size does not match, cannot merge");
        }

        // transform km
        kelvin_mesh.kelvinTransformCentered();
        // kelvin_mesh.plot();

        // add kelvin mesh vertices to this mesh
        for (uint64 i = 0; i < kelvin_mesh.vertices.size(); i++) {
            Vertex* kmv = kelvin_mesh.vertices[i];
            for (uint64_t j = 0; j < boundary_vertices.size(); j++) {
                Vertex* tmb = boundary_vertices[j];
                if (kmv->p == tmb->p) {
                    goto end;
                }
            }
            vertices.push_back(kmv);
            end:;
        }

        // manage neighbors
        for (uint64_t i = 0; i < kelvin_mesh.boundary_vertices.size(); i++) {
            Vertex* kmb = kelvin_mesh.boundary_vertices[i];
            // equivalent mesh boundary vertex
            Vertex* tmb = boundary_vertices[i];
            for (uint8_t j = 0; j < kmb->adjvert_count; j++) {
                // if this adjvert is not a kmb add it to tmb's neighbors
                Vertex* kmb_adj = kmb->adjvert[j];
                bool is_kmb = false;
                for (uint64_t k = 0; k < kelvin_mesh.boundary_vertices.size(); k++) {
                    if (kmb_adj == kelvin_mesh.boundary_vertices[k]) {
                        is_kmb = true;
                        break;
                    }
                }
                if (!is_kmb) {
                    tmb->addAdjacentVertex(kmb_adj);
                    kmb_adj->addAdjacentVertex(tmb);
                }
            }
        }

        /* checks */
        // check kmb in km number
        uint64 kmb_in_km = 0;
        for (uint64_t i = 0; i < kelvin_mesh.boundary_vertices.size(); i++) {
            Vertex* kmb = kelvin_mesh.boundary_vertices[i];
            bool found = false;
            for (uint64_t j = 0; j < kelvin_mesh.vertices.size(); j++) {
                Vertex* kmv = kelvin_mesh.vertices[j];
                if (kmb == kmv) {
                    found = true;
                    kmb_in_km++;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Kelvin mesh boundary vertex not found in kelvin mesh");
            }
        }
        printf("Kelvin mesh boundary vertices in kelvin mesh: %lu\n", kmb_in_km);

        // check if every vertex in the kelvin mesh is in this mesh minus the kelvin boundary vertices
        for (uint64_t i = 0; i < kelvin_mesh.vertices.size(); i++) {
            Vertex* kmv = kelvin_mesh.vertices[i];
            bool is_kmb = false;
            for (uint64_t i = 0; i < kelvin_mesh.boundary_vertices.size(); i++) {
                if (kmv == kelvin_mesh.boundary_vertices[i]) {
                    is_kmb = true;
                    break;
                }
            }
            if (!is_kmb) {
                bool is_in = false;
                for (uint64_t i = 0; i < vertices.size(); i++) {
                    if (kmv == vertices[i]) {
                        is_in = true;
                        break;
                    }
                }
                if (!is_in) {
                    // throw std::runtime_error("Kelvin mesh vertex not in this mesh");
                }
            }
        }

        // check if kelvin boundary vertices are not in this mesh
        for (auto it = kelvin_mesh.boundary_vertices.begin(); it != kelvin_mesh.boundary_vertices.end(); it++) {
            Vertex* kmb = *it;
            for (auto it2 = vertices.begin(); it2 != vertices.end(); it2++) {
                Vertex* tmv = *it2;
                if (kmb == tmv) {
                    throw std::runtime_error("Kelvin boundary vertex found in this mesh");
                }
            }
        }

        // check if every vertex in this mesh has no neighbor in the kelvin mesh boundary
        for (uint64_t i = 0; i < vertices.size(); i++) {
            for (uint8_t j = 0; j < vertices[i]->adjvert_count; j++) {
                Vertex* adj = vertices[i]->adjvert[j];
                for (uint64_t k = 0; k < kelvin_mesh.boundary_vertices.size(); k++) {
                    if (adj == kelvin_mesh.boundary_vertices[k]) {
                        throw std::runtime_error("Vertex in this mesh has neighbor in kelvin mesh boundary");
                    }
                }
            }
        }

        plot();
    }

    void Mesh::enumerateVertices() {
        uint64_t i = 0;
        for (auto it = vertices.begin(); it != vertices.end(); it++) {
            (*it)->id = i;
            i++;
        }
        printf("Enumerated %lu vertices\n", i);
    }

    void Mesh::getFemMatrix(MatCOO &coo) {
        auto start = std::chrono::high_resolution_clock::now();
        for (auto v : vertices) {
            printf("Vertex %lu->", v->id);
            for (uint8_t i = 0; i < v->adjvert_count; i++) {
                printf("%lu,", v->adjvert[i]->id);
            }
            printf("\n");
            uint8_t N = v->adjvert_count;
            cv::Mat S = cv::Mat(N, 5, CV_64F);
            double xi = v->p.x;
            double yi = v->p.y;
            for (uint8_t i = 0; i < N; i++) {
                double xj = v->adjvert[i]->p.x;
                double yj = v->adjvert[i]->p.y;

                double xjmxi = xj - xi;
                double yjmyi = yj - xi;

                S.at<double>(i, 0) = xjmxi;
                S.at<double>(i, 1) = yjmyi;
                S.at<double>(i, 2) = xjmxi * yjmyi;
                S.at<double>(i, 3) = xjmxi * xjmxi * 0.5;
                S.at<double>(i, 4) = yjmyi * yjmyi * 0.5;
            }
            cv::invert(S, S, cv::DECOMP_SVD);
            // print the matrix
            double coeff = 0;
            for (uint8_t i = 0; i < N; i++) {
                // fi * sum c5j - c4j
                coeff += S.at<double>(i, 4) - S.at<double>(i, 3);
            }
            coo.add_elem(v->id, v->id, coeff);
            if (v->id > vertices.size()) {
                printf("Vertex %lu x: %f y: %f\n", v->id, v->p.x, v->p.y);
                throw std::runtime_error("Vertex id out of bounds");
            }
            for (uint8_t i = 0; i < N; i++) {
                // fj * (c4j - c5j)
                coo.add_elem(v->id, v->adjvert[i]->id, S.at<double>(i, 3) - S.at<double>(i, 4));
                if (v->id > vertices.size()) {
                    printf("Vertex %lu x: %f y: %f\n", v->id, v->p.x, v->p.y);
                    throw std::runtime_error("Vertex id out of bounds");
                } else if (v->adjvert[i]->id > vertices.size()) {
                    printf("Vertex %lu x: %f y: %f\n", v->adjvert[i]->id, v->adjvert[i]->p.x, v->adjvert[i]->p.y);
                    throw std::runtime_error("Vertex id out of bounds");
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "FEM matrix construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    }

    void Mesh::getCoefficientVector(CV &b) {
        for (auto v : vertices) {
            uint8_t N = v->adjvert_count;
            double mu_r;
            for (uint8_t i = 0; i < N; i++) {
                mu_r += v->adjmu_r[i];
            }
            mu_r /= N;
            v->mu_r = mu_r;
        }
    }
}