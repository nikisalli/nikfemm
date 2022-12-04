#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <set>

#include "../constants.hpp"
#include "simulation.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/point.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    MagnetostaticSimulation::MagnetostaticSimulation() {

    }

    MagnetostaticSimulation::~MagnetostaticSimulation() {

    }

    void MagnetostaticSimulation::Aplot() {
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL: %s\n");
        }
        SDL_Window* win = SDL_CreateWindow("GAME", // creates a window
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        2000, 2000, 0);
            
        // triggers the program that controls
        // your graphics hardware and sets flags
        // Uint32 render_flags = SDL_RENDERER_ACCELERATED;
    
        // creates a renderer to render our images
        SDL_Renderer* rend = SDL_CreateRenderer(win, -1, 0);

        // clears the window
        SDL_RenderClear(rend);

        // get mesh enclosing rectangle

        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * 2000 - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * 2000 - 0.5 * (max_y + min_y) * y_scale;

        // find max and min
        /*
        float max_A = std::numeric_limits<float>::min();
        float min_A = std::numeric_limits<float>::max();

        for (uint32_t i = 0; i < A.m; i++) {
            if (A[i] > max_A) {
                max_A = A[i];
            }
            if (A[i] < min_A) {
                min_A = A[i];
            }
        }
        */

        std::vector<float> A_sorted(A.val.size());
        for (uint32_t i = 0; i < A.val.size(); i++) {
            A_sorted[i] = A[i];
        }
        std::sort(A_sorted.begin(), A_sorted.end());
        float max_A = A_sorted[0.97 * A_sorted.size()];
        float min_A = A_sorted[0.03 * A_sorted.size()];

        // render

        // clears the window
        SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
        SDL_RenderClear(rend);
        // draw the points
        for (int64_t i = 0; i < mesh.data.numberofpoints; i++) {
            
            // SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
            // // draw a square centered at the point
            // SDL_Rect rect;
            // rect.x = x_offset + v->p.x * x_scale - 2;
            // rect.y = y_offset + v->p.y * y_scale - 2;
            // rect.w = 4;
            // rect.h = 4;
            // SDL_RenderFillRect(rend, &rect);
            // SDL_SetRenderDrawColor(rend, 243, 255, 160, 255);
            // for (uint16_t i = 0; i < v->adjvert_count; i++) {
            //     // SDL_RenderDrawLine(rend, x_offset + v->p.x * x_scale, y_offset + v->p.y * y_scale, x_offset + v->adjvert[i]->p.x * x_scale, y_offset + v->adjvert[i]->p.y * y_scale);
            // }
            
            Point p = mesh.data.pointlist[i];

            if (Point::distance(p, Point(0, 0)) > mesh.radius + mesh.epsilon) {
                continue;
            }

            // verts
            auto points = std::vector<SDL_Vertex>();

            // fill the Vertex<MagnetostaticProp> cell with jet color of the Vertex<MagnetostaticProp>
            SDL_Color c = val2jet(A[i], min_A, max_A);
            // printf("A: %f, c: %d, %d, %d\n", A.val[i], c.r, c.g, c.b);
                            
            // find the triangles that contain the Vertex<MagnetostaticProp> and then
            // for every triangle find the barycenter and add it to the points vector
            for (int32_t j = 0; j < mesh.data.numberoftriangles; j++) {
                if (mesh.data.trianglelist[j].verts[0] == i || mesh.data.trianglelist[j].verts[1] == i || mesh.data.trianglelist[j].verts[2] == i) {
                    SDL_Vertex new_v;
                    Point barycenter = {
                        (mesh.data.pointlist[mesh.data.trianglelist[j].verts[0]].x + mesh.data.pointlist[mesh.data.trianglelist[j].verts[1]].x + mesh.data.pointlist[mesh.data.trianglelist[j].verts[2]].x) / 3,
                        (mesh.data.pointlist[mesh.data.trianglelist[j].verts[0]].y + mesh.data.pointlist[mesh.data.trianglelist[j].verts[1]].y + mesh.data.pointlist[mesh.data.trianglelist[j].verts[2]].y) / 3
                    };
                    new_v.position.x = x_offset + barycenter.x * x_scale;
                    new_v.position.y = y_offset + barycenter.y * y_scale;
                    new_v.color = c;

                    points.push_back(new_v);
                }
            }
            
            // find the center of the points
            SDL_FPoint center = {0, 0};
            for (uint8_t j = 0; j < points.size(); j++) {
                center.x += points[j].position.x;
                center.y += points[j].position.y;
            }
            // printf("count: %d\n", points.size());
            center.x /= points.size();
            center.y /= points.size();
            // printf("center: %d %d\n", center.x, center.y);

            // sort the points by angle
            std::sort(points.begin(), points.end(), [center](SDL_Vertex a, SDL_Vertex b) {
                return atan2(a.position.y - center.y, a.position.x - center.x) < atan2(b.position.y - center.y, b.position.x - center.x);
            });

            // print the points
            
            // for (uint8_t i = 0; i < points.size(); i++) {
            //     printf("%f %f\n", points[i].position.x, points[i].position.y);
            // }
            

            SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
            for (uint8_t j = 0; j < points.size(); j++) {
                SDL_Vertex v1 = points[j];
                SDL_Vertex v2 = points[(j + 1) % points.size()];
                // draw triangle between center and two points
                SDL_Vertex verts[3] = {
                    {center.x, center.y, c.r, c.g, c.b, c.a},
                    {v1.position.x, v1.position.y, c.r, c.g, c.b, c.a},
                    {v2.position.x, v2.position.y, c.r, c.g, c.b, c.a}
                };
                SDL_RenderGeometry(rend, nullptr, verts, 3, nullptr, 0);
            }
        }

        SDL_RenderPresent(rend);

        while(true){
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
            } else if (key.keysym.sym == SDLK_r) {
                x_scale = 1;
                y_scale = 1;
                x_offset = 0;
                y_offset = 0;
            }

            // sleep for 1 second
        }
    }

    void MagnetostaticSimulation::Bplot() {
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL: %s\n");
        }
        SDL_Window* win = SDL_CreateWindow("GAME", // creates a window
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        2000, 2000, 0);

        // triggers the program that controls
        // your graphics hardware and sets flags
        // Uint32 render_flags = SDL_RENDERER_ACCELERATED;
    
        // creates a renderer to render our images
        SDL_Renderer* rend = SDL_CreateRenderer(win, -1, 0);

        // clears the window
        SDL_RenderClear(rend);

        // get mesh enclosing rectangle

        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * 2000 - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * 2000 - 0.5 * (max_y + min_y) * y_scale;

        // find max and min B
        /*
        double max_B = std::numeric_limits<float>::min();
        double min_B = std::numeric_limits<float>::max();

        for (uint32_t i = 0; i < B.size(); i++) {
            double B_mag = B[i].magnitude();
            if (B_mag > max_B) {
                max_B = B_mag;
            }
            if (B_mag < min_B) {
                min_B = B_mag;
            }
        }
        */

        std::vector<double> B_mags;
        for (uint32_t i = 0; i < B.size(); i++) {
            B_mags.push_back(B[i].magnitude());
        }
        std::sort(B_mags.begin(), B_mags.end());
        double max_B = B_mags[0.97 * B_mags.size()];
        double min_B = B_mags[0.03 * B_mags.size()];

        printf("max B: %f\n", max_B);
        printf("min B: %f\n", min_B);

        // max_B = 1e-6;

        // render

        // clears the window
        SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
        SDL_RenderClear(rend);
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Point v1 = mesh.data.pointlist[e.verts[0]];
            Point v2 = mesh.data.pointlist[e.verts[1]];
            Point v3 = mesh.data.pointlist[e.verts[2]];

            if (Point::distance(v1, Point(0, 0)) > mesh.radius + mesh.epsilon ||
                Point::distance(v2, Point(0, 0)) > mesh.radius + mesh.epsilon || 
                Point::distance(v3, Point(0, 0)) > mesh.radius + mesh.epsilon) {
                continue;
            }

            // get the B vectors
            Vector Bv = B[i];
            SDL_Color c = val2jet(Bv.magnitude(), min_B, max_B);
            // get the vertices in window coordinates
            SDL_Vertex verts[3] = {
                {x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset, c.r, c.g, c.b, c.a},
                {x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset, c.r, c.g, c.b, c.a},
                {x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset, c.r, c.g, c.b, c.a}
            };
            // draw the triangle
            SDL_RenderGeometry(rend, nullptr, verts, 3, nullptr, 0);

            // draw mesh edges
            // draw lines from v1 to v2
            SDL_SetRenderDrawColor(rend, 0, 0, 0, 127);
            SDL_RenderDrawLine(rend, x_scale * static_cast<float>(v1.x) + x_offset, 
                                     y_scale * static_cast<float>(v1.y) + y_offset, 
                                     x_scale * static_cast<float>(v2.x) + x_offset, 
                                     y_scale * static_cast<float>(v2.y) + y_offset);
            // draw lines from v2 to v3
            SDL_RenderDrawLine(rend, x_scale * static_cast<float>(v2.x) + x_offset, 
                                     y_scale * static_cast<float>(v2.y) + y_offset, 
                                     x_scale * static_cast<float>(v3.x) + x_offset, 
                                     y_scale * static_cast<float>(v3.y) + y_offset);
            // draw lines from v3 to v1
            SDL_RenderDrawLine(rend, x_scale * static_cast<float>(v3.x) + x_offset, 
                                     y_scale * static_cast<float>(v3.y) + y_offset, 
                                     x_scale * static_cast<float>(v1.x) + x_offset, 
                                     y_scale * static_cast<float>(v1.y) + y_offset);
        }

        // draw the geometry
        // draw the segments
        SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
        for (DrawingSegment s : mesh.drawing.segments) {
            SDL_RenderDrawLine(rend, x_scale * mesh.drawing.points[s.p1].x + x_offset,
                                     y_scale * mesh.drawing.points[s.p1].y + y_offset, 
                                     x_scale * mesh.drawing.points[s.p2].x + x_offset, 
                                     y_scale * mesh.drawing.points[s.p2].y + y_offset);
        }


        // draw the points
        // for (auto v : drawing.points) {
        //     SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
        //     // draw a square centered at the point
        //     SDL_Rect rect;
        //     rect.x = x_offset + v.x * x_scale - 2;
        //     rect.y = y_offset + v.y * y_scale - 2;
        //     rect.w = 4;
        //     rect.h = 4;
        //     SDL_RenderFillRect(rend, &rect);
        // }

        // draw the regions
        // for (DrawingRegion r : drawing.regions) {
        //     SDL_SetRenderDrawColor(rend, 0, 0, 255, 255);
        //     // draw a square centered at the point
        //     SDL_Rect rect;
        //     rect.x = x_offset + r.first.x * x_scale - 4;
        //     rect.y = y_offset + r.first.y * y_scale - 4;
        //     rect.w = 8;
        //     rect.h = 8;
        //     SDL_RenderFillRect(rend, &rect);
        // }

        SDL_RenderPresent(rend);

        while(true){
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
            } else if (key.keysym.sym == SDLK_r) {
                x_scale = 1;
                y_scale = 1;
                x_offset = 0;
                y_offset = 0;
            }

            // sleep for 1 second
        }
    }

    void MagnetostaticSimulation::solve() {
        // get time in milliseconds
        mesh.drawing.addRefiningPoints();

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        auto start1 = std::chrono::high_resolution_clock::now();
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(mesh.drawing.points);
        if (mesh.drawing.points.size() == 0) {
            smallest_circle.radius = 1;
        }
        // translate everything to the origin
        for (uint32_t i = 0; i < mesh.drawing.points.size(); i++) {
            mesh.drawing.points[i] = Point(mesh.drawing.points[i].x - smallest_circle.center.x, mesh.drawing.points[i].y - smallest_circle.center.y);
        }
        for (uint32_t i = 0; i < mesh.drawing.polygons.size(); i++) {
            for (uint32_t j = 0; j < mesh.drawing.polygons[i].points.size(); j++) {
                mesh.drawing.polygons[i].points[j] = Point(mesh.drawing.polygons[i].points[j].x - smallest_circle.center.x, mesh.drawing.polygons[i].points[j].y - smallest_circle.center.y);
            }
        }
        std::vector<DrawingRegion> translated_regions;
        for (DrawingRegion region : mesh.drawing.regions) {
            translated_regions.push_back(DrawingRegion(Point(region.first.x - smallest_circle.center.x, region.first.y - smallest_circle.center.y), region.second));
        }
        mesh.drawing.regions = translated_regions;
        // set simulation offset and boundary radius
        mesh.center = smallest_circle.center;
        mesh.radius = smallest_circle.radius * 2;
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(Point(0, 0), 2 * smallest_circle.radius);
        double circumferential_length = boundary_circle.circumference();
        uint32_t boundary_points = (uint32_t)((circumferential_length / sqrt((MAX_TRIANGLE_AREA * 4) / sqrt(3))) * 2);
        printf("boundary points: %d\n", boundary_points);
        mesh.drawing.drawCircle(boundary_circle, boundary_points);
        // add region near the edge of the circle
        mesh.drawing.drawRegion(Point(boundary_circle.radius * 0.9, 0), {0, {0, 0}, materials::air});
        // add the boundary 
        // mesh.drawing.plot();
        auto start2 = std::chrono::high_resolution_clock::now();
        mesh.mesh();
        auto start3 = std::chrono::high_resolution_clock::now();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        mesh.addKelvinBoundaryConditions(boundary_points);
        mesh.computeEpsilon();
        auto start4 = std::chrono::high_resolution_clock::now();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        printf("the mesh has %u nodes and %u elements\n", mesh.data.numberofpoints, mesh.data.numberoftriangles);
        MatCOO<MagnetostaticNonLinearExpression> coo(mesh.data.numberofpoints);
        CV b(mesh.data.numberofpoints);
        A = CV(mesh.data.numberofpoints);
        B = std::vector<Vector>(mesh.data.numberoftriangles, {0, 0});
        std::vector<float> mu(mesh.data.numberoftriangles, 0);
        std::vector<const MagnetostaticProp*> props(mesh.data.numberoftriangles);

        // fill props
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            props[i] = mesh.drawing.getRegionPtrFromId(mesh.data.triangleattributelist[i]);
        }

        auto start5 = std::chrono::high_resolution_clock::now();

        mesh.getFemSystem(coo, b);
        auto start6 = std::chrono::high_resolution_clock::now();
        mesh.addDirichletBoundaryConditions(coo, b);
        auto start7 = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        printf("coo matrix m: %lu, n: %lu, elems: %lu\n", coo.m, coo.n, coo.elems.size());
        #endif

        MagnetostaticMatCSRSymmetric FemMat(coo);
        // b.write_to_file("b");
        auto start8 = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        // csr.print();
        // coo.plot();
        // x.print();
        // b.print();
        #endif

        // initialize mu
        for (uint32_t i = 0; i < B.size(); i++) {
            mu[i] = props[i]->getMu(0);
        }
        FemMat.updateMat(mu);

        // check if materials are all linear
        bool all_linear = true;
        for (uint32_t i = 0; i < props.size(); i++) {
            if (!props[i]->isLinear()) {
                all_linear = false;
                break;
            }
        }

        if (all_linear) {
            printf("all materials are linear\n");
            preconditionedSSORConjugateGradientSolver(FemMat, b, A, 1.5, 1e-7, 1000);
        } else {
            printf("nonlinear materials detected, starting non linear newton solver\n");
            CV r(b.val.size());  // residual

            double residual = 1e10;
            for (uint32_t i = 0; i < 500; i++) {
                // conjugateGradientSolver(FemMat, b, A, 1e-7, 10000);
                // preconditionedJacobiConjugateGradientSolver(FemMat, b, A, 1e-7, 1000);
                preconditionedSSORConjugateGradientSolver(FemMat, b, A, 1.5, 1e-7, 50);
                // preconditionedIncompleteCholeskyConjugateGradientSolver(FemMat, b, A, 1e-7, 1000);
                // mesh.Aplot(A);
                // mesh.Bplot(B);

                mesh.computeCurl(B, A);
                // mesh.Bplot(B);

                // check if the solution is correct
                FemMat.updateMu(props, mu, B, residual, i);
                FemMat.updateMat(mu);
                CV::mult(r, FemMat, A);
                CV::sub(r, b, r);
                residual = CV::norm(r);
                if (residual < 1e-7) {
                    printf("Converged in %d iterations\n", i);
                    break;
                }
                // printf("%.17g, %.17g; ", K, residual);
                printf("nonlinear iteration %d, residual: %.17g\n", i, residual);
                // print mu
                // printf("%f,", residual);
                fflush(stdout);
            }
        }

        auto start9 = std::chrono::high_resolution_clock::now();

        mesh.computeCurl(B, A);

        printf("%f translate and fix mesh\n", std::chrono::duration_cast<std::chrono::duration<float>>(start2 - start1).count()*1000);
        printf("%f mesh\n", std::chrono::duration_cast<std::chrono::duration<float>>(start3 - start2).count()*1000);
        printf("%f kelvin boundary conditions\n", std::chrono::duration_cast<std::chrono::duration<float>>(start4 - start3).count()*1000);
        printf("%f allocate vectors b x and coo\n", std::chrono::duration_cast<std::chrono::duration<float>>(start5 - start4).count()*1000);
        printf("%f get fem system\n", std::chrono::duration_cast<std::chrono::duration<float>>(start6 - start5).count()*1000);
        printf("%f add dirichlet boundary conditions\n", std::chrono::duration_cast<std::chrono::duration<float>>(start7 - start6).count()*1000);
        printf("%f convert to csr\n", std::chrono::duration_cast<std::chrono::duration<float>>(start8 - start7).count()*1000);
        printf("%f solve\n", std::chrono::duration_cast<std::chrono::duration<float>>(start9 - start8).count()*1000);
        printf("%f total\n", std::chrono::duration_cast<std::chrono::duration<float>>(start9 - start1).count()*1000);

        // auto integrals = mesh.computeForceIntegrals(B);
        // for (auto& i : integrals) {
        //     printf("force integral: %.17gx + %.17gy\n", i.x, i.y);
        // }
        // mesh.muplot(mu);
        // return;
        /*
        mesh.Bplot();

        // report(&out, 1, 1, 0, 0, 0, 0);

        printf("Number of points: %d\nNumber of boundary vertices: %d\n", mesh.vertices.size(), mesh.boundary_vertices.size());
        */

    }
}