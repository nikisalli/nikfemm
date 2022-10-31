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
#include "opencv2/opencv.hpp"

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
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticMesh : Mesh<MagnetostaticProp> {
        MagnetostaticMesh() {
            // default material property
            default_prop = vacuum_prop;
        }

        ~MagnetostaticMesh() {
            
        }

        void Aplot();
        void Bplot();
        void getFemSystem(MatCOO &coo, CV &b);
        // void computeCurl();
    };

    // templated member functions must be defined in the header file
    void MagnetostaticMesh::Aplot() {
        /*
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

        double min_x = -2 * radius;
        double min_y = -2 * radius;
        double max_x = 2 * radius;
        double max_y = 2 * radius;

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

        double max_A = std::numeric_limits<double>::min();
        double min_A = std::numeric_limits<double>::max();

        for (auto v : vertices) {
            if (v->prop.A > max_A) {
                max_A = v->prop.A;
            }
            if (v->prop.A < min_A) {
                min_A = v->prop.A;
            }
        }

        // render

        uint64_t frame = 0;
        while(true){
            // clears the window
            SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
            SDL_RenderClear(rend);
            // draw the points
            uint64_t i = 0;
            for (auto v : vertices) {
                
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
                

                if (geomDistance(v->p, Point(0, 0)) > radius * 2) {
                    continue;
                }

                // verts
                auto points = std::vector<SDL_Vertex>();
                points.reserve(v->adjvert_count);

                // fill the Vertex<MagnetostaticProp> cell with jet color of the Vertex<MagnetostaticProp>
                SDL_Color c = val2jet(v->prop.A, min_A, max_A);
                // printf("A: %f, c: %d, %d, %d\n", v->A, c.r, c.g, c.b);
                                
                // find the triangles that contain the Vertex<MagnetostaticProp> and then
                // for every triangle find the barycenter and add it to the points vector
                for (uint8_t i = 0; i < v->adjvert_count; i++) {
                    for (uint8_t j = 0; j < v->adjvert[i]->adjvert_count; j++) {
                        if (v->adjvert[i]->adjvert[j] == v) {
                            continue;
                        }
                        for (uint8_t k = 0; k < v->adjvert[i]->adjvert[j]->adjvert_count; k++) {
                            if (v->adjvert[i]->adjvert[j]->adjvert[k] == v) {
                                SDL_Vertex new_v;
                                Point barycenter = {
                                    (v->p.x + v->adjvert[i]->p.x + v->adjvert[i]->adjvert[j]->p.x) / 3,
                                    (v->p.y + v->adjvert[i]->p.y + v->adjvert[i]->adjvert[j]->p.y) / 3
                                };
                                new_v.position.x = x_offset + barycenter.x * x_scale;
                                new_v.position.y = y_offset + barycenter.y * y_scale;
                                new_v.color = c;
                                // check if already in points
                                bool already_in = false;
                                for (auto p : points) {
                                    if (p.position.x == new_v.position.x && p.position.y == new_v.position.y) {
                                        already_in = true;
                                        break;
                                    }
                                }
                                if (already_in) {
                                    continue;
                                }
                                points.push_back(new_v);
                            }
                        }
                    }
                }
                
                // find the center of the points
                SDL_FPoint center = {0, 0};
                for (uint8_t i = 0; i < points.size(); i++) {
                    center.x += points[i].position.x;
                    center.y += points[i].position.y;
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
                for (uint8_t i = 0; i < points.size(); i++) {
                    SDL_Vertex v1 = points[i];
                    SDL_Vertex v2 = points[(i + 1) % points.size()];
                    // draw triangle between center and two points
                    SDL_Vertex verts[3] = {
                        {center.x, center.y, c.r, c.g, c.b, c.a},
                        {v1.position.x, v1.position.y, c.r, c.g, c.b, c.a},
                        {v2.position.x, v2.position.y, c.r, c.g, c.b, c.a}
                    };
                    SDL_RenderGeometry(rend, nullptr, verts, 3, nullptr, 0);
                }
                if (i == frame) {
                    // break;
                }
                i++;
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
            } else if (key.keysym.sym == SDLK_r) {
                x_scale = 1;
                y_scale = 1;
                x_offset = 0;
                y_offset = 0;
            } else if (key.keysym.sym == SDLK_SPACE) {
                frame++;
                printf("frame: %d\n", frame);
                Vertex<MagnetostaticProp>* v = vertices[frame];
                SDL_Color c = val2jet(v->prop.A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->prop.A, c.r, c.g, c.b);
            } else if (key.keysym.sym == SDLK_BACKSPACE) {
                frame--;
                printf("frame: %d\n", frame);
                Vertex<MagnetostaticProp>* v = vertices[frame];
                SDL_Color c = val2jet(v->prop.A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->prop.A, c.r, c.g, c.b);
            }

            // sleep for 1 second
        }
        */
    }

    void MagnetostaticMesh::Bplot() {
        /*
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

        double min_x = -2 * radius;
        double min_y = -2 * radius;
        double max_x = 2 * radius;
        double max_y = 2 * radius;

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

        
        double max_A = std::numeric_limits<double>::min();
        double min_A = std::numeric_limits<double>::max();

        for (auto v : vertices) {
            double b_mod = sqrt(v->prop.B.x * v->prop.B.x + v->prop.B.y * v->prop.B.y);
            // printf("bx: %f by: %f b_mod: %f\n", v->B.x, v->B.y, b_mod);
            if (b_mod > max_A) {
                max_A = b_mod;
            }
            if (b_mod < min_A) {
                min_A = b_mod;
            }
        }

        // printf("max_A: %f min_A: %f\n", max_A, min_A);
        

        // render

        uint64_t frame = 0;
        while(true){
            // clears the window
            SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
            SDL_RenderClear(rend);
            // draw the points
            uint64_t i = 0;
            for (auto v : vertices) {

                if (geomDistance(v->p, Point(0, 0)) > radius * 2) {
                    continue;
                }

                // verts
                auto points = std::vector<SDL_Vertex>();
                points.reserve(v->adjvert_count);

                // fill the Vertex<MagnetostaticProp> cell with jet color of the Vertex<MagnetostaticProp>
                double B_mod = sqrt(v->prop.B.x * v->prop.B.x + v->prop.B.y * v->prop.B.y);
                SDL_Color c = val2jet(B_mod, min_A, max_A);
                // printf("A: %f, c: %d, %d, %d\n", v->A, c.r, c.g, c.b);
                                
                // find the triangles that contain the Vertex<MagnetostaticProp> and then
                // for every triangle find the barycenter and add it to the points vector
                for (uint8_t i = 0; i < v->adjvert_count; i++) {
                    for (uint8_t j = 0; j < v->adjvert[i]->adjvert_count; j++) {
                        if (v->adjvert[i]->adjvert[j] == v) {
                            continue;
                        }
                        for (uint8_t k = 0; k < v->adjvert[i]->adjvert[j]->adjvert_count; k++) {
                            if (v->adjvert[i]->adjvert[j]->adjvert[k] == v) {
                                SDL_Vertex new_v;
                                Point barycenter = {
                                    (v->p.x + v->adjvert[i]->p.x + v->adjvert[i]->adjvert[j]->p.x) / 3,
                                    (v->p.y + v->adjvert[i]->p.y + v->adjvert[i]->adjvert[j]->p.y) / 3
                                };
                                new_v.position.x = x_offset + barycenter.x * x_scale;
                                new_v.position.y = y_offset + barycenter.y * y_scale;
                                new_v.color = c;
                                // check if already in points
                                bool already_in = false;
                                for (auto p : points) {
                                    if (p.position.x == new_v.position.x && p.position.y == new_v.position.y) {
                                        already_in = true;
                                        break;
                                    }
                                }
                                if (already_in) {
                                    continue;
                                }
                                points.push_back(new_v);
                            }
                        }
                    }
                }
                
                // find the center of the points
                SDL_FPoint center = {0, 0};
                for (uint8_t i = 0; i < points.size(); i++) {
                    center.x += points[i].position.x;
                    center.y += points[i].position.y;
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

                SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
                for (uint8_t i = 0; i < points.size(); i++) {
                    SDL_Vertex v1 = points[i];
                    SDL_Vertex v2 = points[(i + 1) % points.size()];
                    // draw triangle between center and two points
                    SDL_Vertex verts[3] = {
                        {center.x, center.y, c.r, c.g, c.b, c.a},
                        {v1.position.x, v1.position.y, c.r, c.g, c.b, c.a},
                        {v2.position.x, v2.position.y, c.r, c.g, c.b, c.a}
                    };
                    SDL_RenderGeometry(rend, nullptr, verts, 3, nullptr, 0);
                }
                if (i == frame) {
                    // break;
                }
                i++;
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
            } else if (key.keysym.sym == SDLK_r) {
                x_scale = 1;
                y_scale = 1;
                x_offset = 0;
                y_offset = 0;
            } else if (key.keysym.sym == SDLK_SPACE) {
                frame++;
                printf("frame: %d\n", frame);
                Vertex<MagnetostaticProp>* v = vertices[frame];
                SDL_Color c = val2jet(v->prop.A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->prop.A, c.r, c.g, c.b);
            } else if (key.keysym.sym == SDLK_BACKSPACE) {
                frame--;
                printf("frame: %d\n", frame);
                Vertex<MagnetostaticProp>* v = vertices[frame];
                SDL_Color c = val2jet(v->prop.A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->prop.A, c.r, c.g, c.b);
            } else if (key.keysym.sym == SDLK_f) {
                max_A += 0.1;
            } else if (key.keysym.sym == SDLK_v) {
                max_A -= 0.1;
            } else if (key.keysym.sym == SDLK_g) {
                min_A += 0.1;
            } else if (key.keysym.sym == SDLK_b) {
                min_A -= 0.1;
            }

            // sleep for 1 second
        }
        */
    }

    void MagnetostaticMesh::getFemSystem(MatCOO &coo, CV &b) {
        auto start = std::chrono::high_resolution_clock::now();

        // print points
        for (uint64 i = 0; i < data.numberofpoints; i++) {
            printf("%d %f %f %d\n", i, data.pointlist[i].x, data.pointlist[i].y, 0);
        }

        printf("--------------------\n");

        // print triangles
        for (uint64 i = 0; i < data.numberoftriangles; i++) {
            printf("%d %d %d %d %d\n", i, data.trianglelist[i].verts[0],
                                          data.trianglelist[i].verts[1],
                                          data.trianglelist[i].verts[2],
                                          data.triangleattributelist[i]);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "FEM matrix construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    }

    /*void MagnetostaticMesh::computeCurl() {
        for (auto v : vertices) {
            // find the vertices and their neighbors that are adjacent to v
            std::set<Vertex<MagnetostaticProp>*> adjverts;
            for (uint8_t i = 0; i < v->adjvert_count; i++) {
                Vertex<MagnetostaticProp>* adj_v = v->adjvert[i];
                adjverts.insert(adj_v);
                for (uint8_t j = 0; j < adj_v->adjvert_count; j++) {
                    Vertex<MagnetostaticProp>* adj_adj_v = adj_v->adjvert[j];
                    if (adj_adj_v != v) {
                        adjverts.insert(adj_adj_v);
                    }
                }
            }

            uint8_t N = adjverts.size();
            cv::Mat S = cv::Mat(N, 5, CV_64F);
            cv::Mat b = cv::Mat(N, 1, CV_64F);
            double xi = v->p.x;
            double yi = v->p.y;
            // printf("xi = %f, yi = %f\n", xi, yi);
            uint64_t i = 0;
            for (auto vj : adjverts) {
                double xj = vj->p.x;
                double yj = vj->p.y;

                double xjmxi = xj - xi;
                double yjmyi = yj - yi;
                // printf("xj = %f, yj = %f - {%f, %f, %f}\n", xj, yj, xjmxi, yjmyi, vj->A - v->A);

                S.at<double>(i, 0) = xjmxi;
                S.at<double>(i, 1) = yjmyi;
                S.at<double>(i, 2) = xjmxi * yjmyi;
                S.at<double>(i, 3) = xjmxi * xjmxi * 0.5;
                S.at<double>(i, 4) = yjmyi * yjmyi * 0.5;

                b.at<double>(i, 0) = vj->prop.A - v->prop.A;
                i++;
            }
            // std::cout << S << std::endl;
            // std::cout << b << std::endl;
            cv::invert(S, S, cv::DECOMP_SVD);
            cv::Mat x = S * b;
            // x contains the coefficients of the polynomial in the following order:
            // ax + by + cxy + dx^2 + ey^2
            // the gradient of the polynomial in 0, 0 is a, b
            v->prop.B.x = x.at<double>(1, 0);
            v->prop.B.y = - x.at<double>(0, 0);
            // printf("B = {%f, %f}\n", -v->B.y, v->B.x);
        }
    }*/
}

#endif