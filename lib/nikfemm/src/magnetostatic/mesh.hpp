#ifndef NIK_MAGNETOSTATICMESH_HPP
#define NIK_MAGNETOSTATICMESH_HPP

#include <unordered_set>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <math.h>
#include <set>

#include "SDL2/SDL.h"

#include "../../lib/triangle/triangle.h"
#include "../triangle/util.h"

#include "../constants.hpp"

#include "../utils/utils.hpp"
#include "../drawing/drawing.hpp"
#include "../mesh/mesh.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"
#include "../constants.hpp"
#include "magnetostatic_algebra.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticMesh : Mesh<MagnetostaticProp> {
        MagnetostaticMesh() {
            // default material property
            default_prop = {0, {0, 0}, materials::air};
        }

        ~MagnetostaticMesh() {
            
        }

        void Aplot(CV& A);
        void Bplot(std::vector<Vector>& B);
        void getFemSystem(MatCOO<MagnetostaticNonLinearExpression>& coo, CV& b);
        void addDirichletBoundaryConditions(MatCOO<MagnetostaticNonLinearExpression>& coo, CV& b);
        void computeCurl(std::vector<Vector>& B, CV& A);
    };

    // templated member functions must be defined in the header file
    void MagnetostaticMesh::Aplot(CV &A) {
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL: %s\n");
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

        float min_x = -2 * radius;
        float min_y = -2 * radius;
        float max_x = 2 * radius;
        float max_y = 2 * radius;

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

        // render

        uint32_t frame = 0;
        // clears the window
        SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
        SDL_RenderClear(rend);
        // draw the points
        uint32_t i = 0;
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            
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
            
            Point p = data.pointlist[i];

            if (geomDistance(p, Point(0, 0)) > radius * 2) {
                continue;
            }

            // verts
            auto points = std::vector<SDL_Vertex>();

            // fill the Vertex<MagnetostaticProp> cell with jet color of the Vertex<MagnetostaticProp>
            SDL_Color c = val2jet(A[i], min_A, max_A);
            // printf("A: %f, c: %d, %d, %d\n", A.val[i], c.r, c.g, c.b);
                            
            // find the triangles that contain the Vertex<MagnetostaticProp> and then
            // for every triangle find the barycenter and add it to the points vector
            for (uint32_t j = 0; j < data.numberoftriangles; j++) {
                if (data.trianglelist[j].verts[0] == i || data.trianglelist[j].verts[1] == i || data.trianglelist[j].verts[2] == i) {
                    SDL_Vertex new_v;
                    Point barycenter = {
                        (data.pointlist[data.trianglelist[j].verts[0]].x + data.pointlist[data.trianglelist[j].verts[1]].x + data.pointlist[data.trianglelist[j].verts[2]].x) / 3,
                        (data.pointlist[data.trianglelist[j].verts[0]].y + data.pointlist[data.trianglelist[j].verts[1]].y + data.pointlist[data.trianglelist[j].verts[2]].y) / 3
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

    void MagnetostaticMesh::Bplot(std::vector<Vector>& B) {
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL: %s\n");
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

        float min_x = -2 * radius;
        float min_y = -2 * radius;
        float max_x = 2 * radius;
        float max_y = 2 * radius;

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

        printf("max B: %f\n", max_B);
        printf("min B: %f\n", min_B);

        // max_B = 1e-6;

        // render

        uint32_t frame = 0;
        // clears the window
        SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
        SDL_RenderClear(rend);
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            // get the triangle
            Elem e = data.trianglelist[i];
            // get the vertices
            Point v1 = data.pointlist[e.verts[0]];
            Point v2 = data.pointlist[e.verts[1]];
            Point v3 = data.pointlist[e.verts[2]];

            if (geomDistance(v1, Point(0, 0)) > 2 * radius || geomDistance(v2, Point(0, 0)) > 2 * radius || geomDistance(v3, Point(0, 0)) > 2 * radius) {
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
        }

        // draw the geometry
        // draw the segments
        SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
        for (DrawingSegment s : drawing.segments) {
            SDL_RenderDrawLine(rend, x_scale * drawing.points[s.p1].x + x_offset, y_scale * drawing.points[s.p1].y + y_offset, x_scale * drawing.points[s.p2].x + x_offset, y_scale * drawing.points[s.p2].y + y_offset);
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

    void MagnetostaticMesh::getFemSystem(MatCOO<MagnetostaticNonLinearExpression>&coo, CV &b) {
        auto start = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            printf("%d %d\n", i, (int)drawing.getRegionFromId(data.triangleattributelist[i]).mu);
        }
        #endif

        auto adjelems_ids = new uint32_t[data.numberofpoints][18];
        auto adjelems_props = new const MagnetostaticProp*[data.numberofpoints][18];
        auto adjelems_count = new uint8_t[data.numberofpoints]();  // initialize to 0

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = data.trianglelist[i].verts[j];
                adjelems_ids[myid][adjelems_count[myid]] = i;
                adjelems_props[myid][adjelems_count[myid]++] = drawing.getRegionPtrFromId(data.triangleattributelist[i]);
            }
        }
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            for (uint8_t j = 0; j < adjelems_count[i]; j++) {
                uint32_t v1, v2, v3;
                v1 = i;
                Elem myelem = data.trianglelist[adjelems_ids[i][j]];
                if (i == data.trianglelist[adjelems_ids[i][j]].verts[0]) {
                    v2 = myelem.verts[1];
                    v3 = myelem.verts[2];
                } else if (i == myelem.verts[1]) {
                    v2 = myelem.verts[2];
                    v3 = myelem.verts[0];
                } else if (i == myelem.verts[2]) {
                    v2 = myelem.verts[0];
                    v3 = myelem.verts[1];
                } else {
                    nexit("error: vertex not found in element");
                }

                double oriented_area = Point::double_oriented_area(data.pointlist[v1], data.pointlist[v2], data.pointlist[v3]);
                
                if (oriented_area < 0) {
                    std::swap(v2, v3);
                }

                double area = Point::double_oriented_area(data.pointlist[v1], data.pointlist[v2], data.pointlist[v3]);

                double b1 = (data.pointlist[v2].y - data.pointlist[v3].y) / area;
                double c1 = (data.pointlist[v3].x - data.pointlist[v2].x) / area;
                // coo.add_elem(i, v1, (area * (b1 * b1 + c1 * c1)) / (2 * adjelems_props[i][j].mu));
                // check if key exists
                uint64_t key = MatCOO<int>::get_key(i, v1);
                if (coo.elems.find(key) == coo.elems.end()) {
                    coo.elems.emplace(key, MagnetostaticNonLinearExpression());
                }
                coo.elems[key].terms.push_back(
                    {
                        (area * (b1 * b1 + c1 * c1) * 0.5),
                        adjelems_ids[i][j],
                        false
                    }
                );
                if (v2 <= i) {
                    double b2 = (data.pointlist[v3].y - data.pointlist[v1].y) / area;
                    double c2 = (data.pointlist[v1].x - data.pointlist[v3].x) / area;
                    // coo.add_elem(i, v2, (area * (b2 * b1 + c2 * c1)) / (2 * adjelems_props[i][j].mu));
                    key = MatCOO<int>::get_key(i, v2);
                    if (coo.elems.find(key) == coo.elems.end()) {
                        coo.elems.emplace(key, MagnetostaticNonLinearExpression());
                    }
                    coo.elems[key].terms.push_back(
                        {
                            (area * (b2 * b1 + c2 * c1) * 0.5),
                            adjelems_ids[i][j],
                            false
                        }
                    );
                }
                if (v3 <= i) {
                    double b3 = (data.pointlist[v1].y - data.pointlist[v2].y) / area;
                    double c3 = (data.pointlist[v2].x - data.pointlist[v1].x) / area;
                    // coo.add_elem(i, v3, (area * (b3 * b1 + c3 * c1)) / (2 * adjelems_props[i][j].mu));
                    key = MatCOO<int>::get_key(i, v3);
                    if (coo.elems.find(key) == coo.elems.end()) {
                        coo.elems.emplace(key, MagnetostaticNonLinearExpression());
                    }
                    coo.elems[MatCOO<int>::get_key(i, v3)].terms.push_back(
                        {
                            (area * (b3 * b1 + c3 * c1) * 0.5),
                            adjelems_ids[i][j],
                            false
                        }
                    );
                }

                // set the b vector
                b.add_elem(i, (area * adjelems_props[i][j]->J) / 6);
            }
        }

        // iterate over upper triangular matrix and copy to lower triangular matrix
        // for (uint32_t i = 0; i < coo.m; i++) {
        //     for (uint32_t j = i; j < coo.n; j++) {
        //         if (coo.get_elem(i, j) != 0) {
        //             coo.add_elem(j, i, coo.get_elem(i, j));
        //         }
        //     }
        // }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "FEM matrix construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    }

    void MagnetostaticMesh::addDirichletBoundaryConditions(MatCOO<MagnetostaticNonLinearExpression> &coo, CV &b) {
        // find three furthest points from the center
        uint32_t p1, p2, p3;
        double d1, d2, d3;
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            double dist = geomDistance(data.pointlist[i], Point(0, 0));
            if (dist > d1) {
                d3 = d2;
                p3 = p2;
                d2 = d1;
                p2 = p1;
                d1 = dist;
                p1 = i;
            } else if (dist > d2) {
                d3 = d2;
                p3 = p2;
                d2 = dist;
                p2 = i;
            } else if (dist > d3) {
                d3 = dist;
                p3 = i;
            }
        }

        for (auto elem : coo.elems) {
            uint32_t m = elem.first >> 32;
            uint32_t n = elem.first & 0xFFFFFFFF;
            if (m == p1 || m == p2 || m == p3) {
                coo.elems[MatCOO<int>::get_key(m, n)].terms.clear();
                coo.elems[MatCOO<int>::get_key(m, n)].terms.push_back({0, 0, true});
            }
            if (n == p1 || n == p2 || n == p3) {
                coo.elems[MatCOO<int>::get_key(m, n)].terms.clear();
                coo.elems[MatCOO<int>::get_key(m, n)].terms.push_back({0, 0, true});
            }
        }

        coo.elems[MatCOO<int>::get_key(p1, p1)].terms.clear();
        coo.elems[MatCOO<int>::get_key(p1, p1)].terms.push_back({ 1, 0, true});
        coo.elems[MatCOO<int>::get_key(p2, p2)].terms.clear();
        coo.elems[MatCOO<int>::get_key(p2, p2)].terms.push_back({ 1, 0, true});
        coo.elems[MatCOO<int>::get_key(p3, p3)].terms.clear();
        coo.elems[MatCOO<int>::get_key(p3, p3)].terms.push_back({ 1, 0, true});
        // since we are setting the potential to zero at the three points, we do not need to subtract anything from the b vector
    }

    void MagnetostaticMesh::computeCurl(std::vector<Vector>& B, CV &A) {
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            Elem myelem = data.trianglelist[i];
            double x1 = data.pointlist[myelem.verts[0]].x;
            double y1 = data.pointlist[myelem.verts[0]].y;
            double z1 = A[myelem.verts[0]];
            double x2 = data.pointlist[myelem.verts[1]].x;
            double y2 = data.pointlist[myelem.verts[1]].y;
            double z2 = A[myelem.verts[1]];
            double x3 = data.pointlist[myelem.verts[2]].x;
            double y3 = data.pointlist[myelem.verts[2]].y;
            double z3 = A[myelem.verts[2]];

            // fit z = a + bx + cy to the three points
            double a1 = x2 - x1;
            double b1 = y2 - y1;
            double c1 = z2 - z1;
            double a2 = x3 - x1;
            double b2 = y3 - y1;
            double c2 = z3 - z1;

            double a = b1 * c2 - b2 * c1;
            double b = a2 * c1 - a1 * c2;
            double c = a1 * b2 - b1 * a2;

            double dx = -a / c;
            double dy = -b / c;
            
            // printf("dx = %.17g dy = %.17g for elem (%.1f, %.1f, %.17g), (%.1f, %.1f, %.17g), (%.1f, %.1f, %.17g)\n", dx, dy, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            B[i] = Vector(-dy, dx);
        }
    }
}

#endif