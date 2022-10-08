#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <math.h>
#include <set>

#include "SDL2/SDL.h"
#include "SDL2/SDL_ttf.h"
#include "opencv2/opencv.hpp"

#include "../../lib/triangle/triangle.h"
#include "../triangle/util.h"

#include <constants.hpp>

#include "mesh.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    Mesh::Mesh() {
        
    }

    Mesh::~Mesh() {

    }

    void Mesh::plot() {
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

    void Mesh::Aplot() {
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
        /*
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
        */

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
            if (v->A > max_A) {
                max_A = v->A;
            }
            if (v->A < min_A) {
                min_A = v->A;
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
                /*
                SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
                // draw a square centered at the point
                SDL_Rect rect;
                rect.x = x_offset + v->p.x * x_scale - 2;
                rect.y = y_offset + v->p.y * y_scale - 2;
                rect.w = 4;
                rect.h = 4;
                SDL_RenderFillRect(rend, &rect);

                SDL_SetRenderDrawColor(rend, 243, 255, 160, 255);
                for (uint16_t i = 0; i < v->adjvert_count; i++) {
                    // SDL_RenderDrawLine(rend, x_offset + v->p.x * x_scale, y_offset + v->p.y * y_scale, x_offset + v->adjvert[i]->p.x * x_scale, y_offset + v->adjvert[i]->p.y * y_scale);
                }
                */

                if (geomDistance(v->p, Point(0, 0)) > radius * 2) {
                    continue;
                }

                // verts
                auto points = std::vector<SDL_Vertex>();
                points.reserve(v->adjvert_count);

                // fill the vertex cell with jet color of the vertex
                SDL_Color c = val2jet(v->A, min_A, max_A);
                // printf("A: %f, c: %d, %d, %d\n", v->A, c.r, c.g, c.b);
                                
                // find the triangles that contain the vertex and then
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
                /*
                for (uint8_t i = 0; i < points.size(); i++) {
                    printf("%f %f\n", points[i].position.x, points[i].position.y);
                }
                */

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
                Vertex* v = vertices[frame];
                SDL_Color c = val2jet(v->A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->A, c.r, c.g, c.b);
            } else if (key.keysym.sym == SDLK_BACKSPACE) {
                frame--;
                printf("frame: %d\n", frame);
                Vertex* v = vertices[frame];
                SDL_Color c = val2jet(v->A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->A, c.r, c.g, c.b);
            }

            // sleep for 1 second
        }
    }

    void Mesh::Bplot() {
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
        /*
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::min();
        double max_y = std::numeric_limits<double>::max();
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
        */

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
            double b_mod = sqrt(v->B.x * v->B.x + v->B.y * v->B.y);
            // printf("bx: %f by: %f b_mod: %f\n", v->B.x, v->B.y, b_mod);
            if (b_mod > max_A) {
                max_A = b_mod;
            }
            if (b_mod < min_A) {
                min_A = b_mod;
            }
        }

        // printf("max_A: %f min_A: %f\n", max_A, min_A);
        

        /*
        double average_A = 0;
        for (auto v : vertices) {
            double b_mod = sqrt(v->B.x * v->B.x + v->B.y * v->B.y);
            average_A += b_mod;
        }

        average_A /= vertices.size();
        double max_A = 2 * average_A;
        double min_A = 0;
        */

        // render

        uint64_t frame = 0;
        while(true){
            // clears the window
            SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
            SDL_RenderClear(rend);
            // draw the points
            uint64_t i = 0;
            for (auto v : vertices) {
                /*
                SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
                // draw a square centered at the point
                SDL_Rect rect;
                rect.x = x_offset + v->p.x * x_scale - 2;
                rect.y = y_offset + v->p.y * y_scale - 2;
                rect.w = 4;
                rect.h = 4;
                SDL_RenderFillRect(rend, &rect);

                SDL_SetRenderDrawColor(rend, 243, 255, 160, 255);
                for (uint16_t i = 0; i < v->adjvert_count; i++) {
                    // SDL_RenderDrawLine(rend, x_offset + v->p.x * x_scale, y_offset + v->p.y * y_scale, x_offset + v->adjvert[i]->p.x * x_scale, y_offset + v->adjvert[i]->p.y * y_scale);
                }
                */

                if (geomDistance(v->p, Point(0, 0)) > radius * 2) {
                    continue;
                }

                // verts
                auto points = std::vector<SDL_Vertex>();
                points.reserve(v->adjvert_count);

                // fill the vertex cell with jet color of the vertex
                double B_mod = sqrt(v->B.x * v->B.x + v->B.y * v->B.y);
                SDL_Color c = val2jet(B_mod, min_A, max_A);
                // printf("A: %f, c: %d, %d, %d\n", v->A, c.r, c.g, c.b);
                                
                // find the triangles that contain the vertex and then
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
                /*
                for (uint8_t i = 0; i < points.size(); i++) {
                    printf("%f %f\n", points[i].position.x, points[i].position.y);
                }
                */

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
                Vertex* v = vertices[frame];
                SDL_Color c = val2jet(v->A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->A, c.r, c.g, c.b);
            } else if (key.keysym.sym == SDLK_BACKSPACE) {
                frame--;
                printf("frame: %d\n", frame);
                Vertex* v = vertices[frame];
                SDL_Color c = val2jet(v->A, min_A, max_A);
                printf("A: %f %d %d %d\n", v->A, c.r, c.g, c.b);
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
            // printf("adding region %f %f %lu\n", r.first.x, r.first.y, r.second);
            in.regionlist[4 * i] = r.first.x;
            in.regionlist[4 * i + 1] = r.first.y;
            in.regionlist[4 * i + 2] = r.second;
            in.regionlist[4 * i + 3] = 0;
            i++;
        }
        // printf("----------------------\n");

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

            // add muj
            // v1->addAdjacentMu(out.triangleattributelist[i]);
            // v2->addAdjacentMu(out.triangleattributelist[i]);
            // v3->addAdjacentMu(out.triangleattributelist[i]);

            double reg_val = drawing.getRegionFromId(out.triangleattributelist[i]);
            v1->addAdjacentMu(reg_val);
            v2->addAdjacentMu(reg_val);
            v3->addAdjacentMu(reg_val);
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
        const double max_radius_coeff = 100;
        double R_squared = radius * radius;
        double max_x = radius / max_radius_coeff;
        // printf("R^2 = %f\n", R_squared);

        for (auto it = vertices.begin(); it != vertices.end(); it++) {
            Vertex* v = *it;
            double mag_squared = v->p.x * v->p.x + v->p.y * v->p.y;
            double scale = R_squared / mag_squared;
            // printf("v = (%f, %f) -> (%f, %f), mag = %f, scale = %f\n", v->p.x, v->p.y, v->p.x * scale, v->p.y * scale, mag_squared, scale);
            double dist = geomDistance(v->p, center);
            if (dist < max_x) {
            // if (false) {
#ifdef DEBUG_PRINT
                printf("kelvin transform too large\n");
#endif
                v->p = v->p * ((radius * max_radius_coeff) / dist);
            } else {
                v->p = v->p * scale;
            }
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
            nexit("Kelvin mesh boundary vertices size does not match, cannot merge");
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
        /*
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
                nexit("Kelvin mesh boundary vertex not found in kelvin mesh");
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
                    nexit("Kelvin mesh vertex not in this mesh");
                }
            }
        }

        // check if kelvin boundary vertices are not in this mesh
        for (auto it = kelvin_mesh.boundary_vertices.begin(); it != kelvin_mesh.boundary_vertices.end(); it++) {
            Vertex* kmb = *it;
            for (auto it2 = vertices.begin(); it2 != vertices.end(); it2++) {
                Vertex* tmv = *it2;
                if (kmb == tmv) {
                    nexit("Kelvin boundary vertex found in this mesh");
                }
            }
        }

        // check if every vertex in this mesh has no neighbor in the kelvin mesh boundary
        for (uint64_t i = 0; i < vertices.size(); i++) {
            for (uint8_t j = 0; j < vertices[i]->adjvert_count; j++) {
                Vertex* adj = vertices[i]->adjvert[j];
                for (uint64_t k = 0; k < kelvin_mesh.boundary_vertices.size(); k++) {
                    if (adj == kelvin_mesh.boundary_vertices[k]) {
                        nexit("Vertex in this mesh has neighbor in kelvin mesh boundary");
                    }
                }
            }
        }*/

        // plot();
    }

    void Mesh::addDirichletBoundaryConditions(MatCOO &coo, CV &b) {
        // find the three furthest vertices from the center
        double dist1, dist2, dist3;
        uint64_t index1, index2, index3;

        for (uint64_t i = 0; i < vertices.size(); i++) {
            Vertex* v = vertices[i];
            double dist = geomDistance(v->p, center);
            if (dist > dist1) {
                dist3 = dist2;
                index3 = index2;
                dist2 = dist1;
                index2 = index1;
                dist1 = dist;
                index1 = v->id;
            } else if (dist > dist2 && dist != dist1) {
                dist3 = dist2;
                index3 = index2;
                dist2 = dist;
                index2 = v->id;
            } else if (dist > dist3 && dist != dist2) {
                dist3 = dist;
                index3 = v->id;
            }
        }

        // all the coefficients on the inedxN-th column of the matrix A are set to 0 for dirichlet boundary conditions
        for (uint64_t i = 0; i < coo.elems.size(); i++) {
            ElemCOO e = coo.elems[i];
            if (e.n == index1 || e.n == index2 || e.n == index3 || e.m == index1 || e.m == index2 || e.m == index3) {
                if (e.n == e.m) {
                    e.val = 1;
#ifdef DEBUG_PRINT
                    printf("Dirichlet boundary condition on FEM matrix element %lu %lu (element is on the diagonal)\n", e.m, e.n);
#endif
                } else {
                    e.val = 0;
#ifdef DEBUG_PRINT
                    printf("Dirichlet boundary condition on FEM matrix element %lu %lu\n", e.m, e.n);
#endif
                }

                /*
                // in this case the dirichlet boundary condition sets the magnetic vector potential A to zero on the boundary so it is
                // not necessary to subtract the coefficient of the corresponding value multiplied by the dirichlet magnetic vector potential
                // from the b vector
                if (e.m == index1 || e.m == index2 || e.m == index3) {
                    b[e.m] -= e.val * 0;
                }
                */
            }
            coo.elems[i] = e;
        }
    }

    void Mesh::enumerateVertices() {
        uint64_t i = 0;
        for (auto it = vertices.begin(); it != vertices.end(); it++) {
            (*it)->id = i;
            i++;
        }
#ifdef DEBUG_PRINT
        printf("Enumerated %lu vertices\n", i);
#endif
    }

    void Mesh::getFemMatrix(MatCOO &coo) {
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t num = 0;
        for (auto v : vertices) {
            double sum = 0;
            for (uint8_t i = 0; i < v->adjvert_count; i++) {
                Vertex* adj_v = v->adjvert[i];
                // find the two opposite vertices relative to the edge between v and adj_v
                // to do this we find the two vertices that are not v or adj_v and are adjacent to both v and adj_v
                Vertex* opp_v1 = nullptr;
                Vertex* opp_v2 = nullptr;
                for (uint8_t j = 0; j < v->adjvert_count; j++) {
                    Vertex* v_adj = v->adjvert[j];
                    if (v_adj != adj_v) {
                        for (uint8_t k = 0; k < adj_v->adjvert_count; k++) {
                            Vertex* adj_v_adj = adj_v->adjvert[k];
                            if (adj_v_adj != v) {
                                if (v_adj == adj_v_adj) {
                                    if (opp_v1 == nullptr) {
                                        opp_v1 = v_adj;
                                    } else if (opp_v2 == nullptr) {
                                        opp_v2 = v_adj;
                                    } else {
                                        nexit("More than two opposite vertices found");
                                    }
                                }
                            }
                        }
                    }
                }
                if (opp_v1 == nullptr || opp_v2 == nullptr) {
                    nexit("Less than two opposite vertices found");
                }
                // get angle v-opp_v1-adj_v
                double angle1 = geomAngle(v->p, opp_v1->p, adj_v->p);
                // get angle v-opp_v2-adj_v
                double angle2 = geomAngle(v->p, opp_v2->p, adj_v->p);
                // double w = 0.5 * (1 / tan(angle1) + 1 / tan(angle2));
                double w = 0.5 * (cos(angle1) / sin(angle1) + cos(angle2) / sin(angle2));
#ifdef DEBUG_PRINT
                printf ("elem num %lu\n", num);
                printf ("angle1 = %f\n", angle1);
                printf ("angle2 = %f\n", angle2);
                printf ("w = %f\n", w);
#endif
                // add the coefficient to the matrix
                coo.add_elem(v->id, adj_v->id, w);
                sum += w;
            }
#ifdef DEBUG_PRINT
            printf("sum = %f\n", sum);
            printf("------------------------------------------------------------\n");
#endif
            // add the coefficient to the matrix
            coo.add_elem(v->id, v->id, -sum);
            num++;
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "FEM matrix construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    }

    void Mesh::getCoefficientVector(CV &b) {
        for (auto v : vertices) {
            uint8_t N = v->adjmuj_count;
            double mu_r = 0;
            for (uint8_t i = 0; i < N; i++) {
                mu_r += v->adjmuj[i];
#ifdef DEBUG_PRINT
                printf("%f ", v->adjmuj[i]);
#endif
            }
            mu_r /= N;
#ifdef DEBUG_PRINT
            printf("\nmu_r = %f\n", mu_r);
#endif
            v->muj = mu_r;
            // find area of the cell centered at v
            double area = v->cellArea();
            b[v->id] = mu_r * area;
        }
    }

    void Mesh::setField(CV &x) {
        for (auto v : vertices) {
            v->A = x[v->id];
        }
    }

    void Mesh::computeCurl() {
        for (auto v : vertices) {
            // find the vertices and their neighbors that are adjacent to v
            std::set<Vertex*> adjverts;
            for (uint8_t i = 0; i < v->adjvert_count; i++) {
                Vertex* adj_v = v->adjvert[i];
                adjverts.insert(adj_v);
                for (uint8_t j = 0; j < adj_v->adjvert_count; j++) {
                    Vertex* adj_adj_v = adj_v->adjvert[j];
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

                b.at<double>(i, 0) = vj->A - v->A;
                i++;
            }
            // std::cout << S << std::endl;
            // std::cout << b << std::endl;
            cv::invert(S, S, cv::DECOMP_SVD);
            cv::Mat x = S * b;
            // x contains the coefficients of the polynomial in the following order:
            // ax + by + cxy + dx^2 + ey^2
            // the gradient of the polynomial in 0, 0 is a, b
            v->B.x = x.at<double>(1, 0);
            v->B.y = - x.at<double>(0, 0);
            // printf("B = {%f, %f}\n", -v->B.y, v->B.x);
        }
    }
}