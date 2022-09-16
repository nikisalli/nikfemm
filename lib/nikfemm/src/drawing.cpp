#include <unordered_set>
#include <set>
#include <chrono>
#include <array>

#include "SDL2/SDL.h"

#include <constants.hpp>

#include "drawing.hpp"

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

    void Drawing::plot() {
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

        std::unordered_set<Point> vertices;

        for (Segment s : segments) {
            vertices.insert(s.p1);
            vertices.insert(s.p2);
        }

        // get mesh enclosing rectangle
        double min_x = 1000000;
        double min_y = 1000000;
        double max_x = -1000000;
        double max_y = -1000000;
        for (auto v : vertices) {
            if (v.x < min_x) {
                min_x = v.x;
            }
            if (v.y < min_y) {
                min_y = v.y;
            }
            if (v.x > max_x) {
                max_x = v.x;
            }
            if (v.y > max_y) {
                max_y = v.y;
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
            // draw the segments
            for (Segment s : segments) {
                SDL_RenderDrawLine(rend, x_offset + s.p1.x * x_scale, y_offset + s.p1.y * y_scale, x_offset + s.p2.x * x_scale, y_offset + s.p2.y * y_scale);
            }
    

            // draw the points
            for (auto v : vertices) {
                SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
                // draw a square centered at the point
                SDL_Rect rect;
                rect.x = x_offset + v.x * x_scale - 2;
                rect.y = y_offset + v.y * y_scale - 2;
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