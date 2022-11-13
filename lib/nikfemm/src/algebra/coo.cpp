#include <math.h>

#include "assert.h"

#include "SDL2/SDL.h"

#include "coo.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    MatCOO::MatCOO(uint32_t m) {
        this->m = m;
    }

    MatCOO::~MatCOO() {
    }

    double MatCOO::operator()(uint32_t _m, uint32_t _n) const {
        // swap indices if _m > _n
        if (_m > _n) {
            std::swap(_m, _n);
        }
        uint64_t key = (uint64_t)_m << 32 | (uint64_t)_n;
        if (elems.find(key) != elems.end()) {
            return elems.at(key);
        } else {
            return 0.0;
        }
    }

    void MatCOO::set_elem(uint32_t _m, uint32_t _n, double val) {
        // assert(_m < m);
        // assert(_n < n);
        // assert(_m <= n);
        // if _m > _n, swap indices
        if (_m > _n) {
            std::swap(_m, _n);
        }
        elems[(uint64_t)_m << 32 | (uint64_t)_n] = val;
    }

    void MatCOO::add_elem(uint32_t _m, uint32_t _n, double val) {
        // is element in matrix?
        // assert(_m < m);
        // assert(_n < n);
        // assert(_m <= n);
        // if _m > _n, swap indices
        if (_m > _n) {
            std::swap(_m, _n);
        }
        if (elems.find((uint64_t)_m << 32 | (uint64_t)_n) != elems.end()) {
            elems[(uint64_t)_m << 32 | (uint64_t)_n] += val;
        } else {
            elems[(uint64_t)_m << 32 | (uint64_t)_n] = val;
        }
    }

    double MatCOO::get_elem(uint32_t _m, uint32_t _n) {
        // assert(_m < m);
        // assert(_n < n);
        // assert(_m <= n);
        // if _m > _n, swap indices
        if (_m > _n) {
            std::swap(_m, _n);
        }
        if (elems.find((uint64_t)_m << 32 | (uint64_t)_n) != elems.end()) {
            return elems[(uint64_t)_m << 32 | (uint64_t)_n];
        } else {
            return 0.0;
        }
    }

    void MatCOO::print() {
        uint32_t prev_m = 0;
        uint32_t prev_n = 0;
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
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
                    // plot the symmetrical element with respect to the diagonal
                    rect.x = x_offset + _m * x_scale;
                    rect.y = y_offset + _n * y_scale;
                    SDL_RenderFillRect(rend, &rect);
                }
            }

            // draw square around the matrix
            SDL_Rect rect;
            rect.x = x_offset;
            rect.y = y_offset;
            rect.w = m * x_scale;
            rect.h = m * y_scale;
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