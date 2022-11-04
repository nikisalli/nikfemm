#include <math.h>

#include "SDL2/SDL.h"

#include "coo.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    MatCOO::MatCOO() {
        m = 0;
        n = 0;
    }

    MatCOO::~MatCOO() {
    }

    void MatCOO::set_elem(uint32_t _m, uint32_t _n, double val) {
        if (_m + 1 > m) m = _m + 1;
        if (_n + 1 > n) n = _n + 1;
        elems[(uint64_t)_m << 32 | (uint64_t)_n] = val;
    }

    void MatCOO::add_elem(uint32_t _m, uint32_t _n, double val) {
        if (_m + 1 > m) m = _m + 1;
        if (_n + 1 > n) n = _n + 1;
        // is element in matrix?
        if (elems.find((uint64_t)_m << 32 | (uint64_t)_n) != elems.end()) {
            elems[(uint64_t)_m << 32 | (uint64_t)_n] += val;
        } else {
            elems[(uint64_t)_m << 32 | (uint64_t)_n] = val;
        }
    }

    double MatCOO::get_elem(uint32_t _m, uint32_t _n) {
        if (elems.find((uint64_t)_m << 32 | (uint64_t)_n) != elems.end()) {
            return elems[(uint64_t)_m << 32 | (uint64_t)_n];
        } else {
            return 0.0;
        }
    }

    void MatCOO::plot() {
        // draw the mesh
        // returns zero on success else non-zero
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL");
            exit(1);
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

        double x_scale = 1;
        double y_scale = 1;
        double x_offset = 0;
        double y_offset = 0;

        // draw the matrix as a grid of squares where if the value is non-zero, the square is filled. make it zoomable and pannable
        while(true){
            SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
            SDL_RenderClear(rend);
            SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
            for (auto elem : elems) {
                uint32_t _m = elem.first >> 32;
                uint32_t _n = elem.first & 0xFFFFFFFF;
                double val = elem.second;
                if (val != 0) {
                    SDL_Rect rect;
                    rect.x = x_offset + _n * x_scale;
                    rect.y = y_offset + _m * y_scale;
                    rect.w = x_scale;
                    rect.h = y_scale;
                    SDL_RenderFillRect(rend, &rect);
                }
            }

            // draw square around the matrix
            SDL_Rect rect;
            rect.x = x_offset;
            rect.y = y_offset;
            rect.w = m * x_scale;
            rect.h = n * y_scale;
            SDL_RenderDrawRect(rend, &rect);

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
        }
    }
}