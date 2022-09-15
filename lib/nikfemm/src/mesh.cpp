#include "SDL2/SDL.h"

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
}