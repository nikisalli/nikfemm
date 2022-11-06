#ifndef NIK_DRAWING_HPP
#define NIK_DRAWING_HPP

#include <unordered_set>
#include <utility>
#include <vector>
#include <map>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <set>
#include <chrono>
#include <array>

#include "SDL2/SDL.h"

#include "../constants.hpp"

#include "../utils/utils.hpp"
#include "../geometry/point.hpp"
#include "../geometry/segment.hpp"
#include "drawing_segment.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"

namespace nikfemm {
    typedef std::pair<Point, uint32_t> DrawingRegion;

    template<typename Prop>
    struct Drawing {
        std::map<Prop, uint32_t> region_map;
        std::vector<DrawingRegion> regions;
        std::vector<Point> points;
        std::unordered_set<DrawingSegment> segments;

        Drawing();
        ~Drawing();

        /* drawing */
        public:
            void drawRectangle(Point p1, Point p2);
            void drawRectangle(Point p1, double width, double height);
            void drawCircle(Point p, double radius, uint32_t n_segments);
            void drawCircle(Circle c, uint32_t n_segments);
            void drawPolygon(Point* points, uint32_t n_points);
            void drawPolyLine(Point* points, uint32_t n_points);  // same as drawPolygon, but doesn't close the figure
            void drawRegion(Point p, Prop val);
            void drawSegment(Point p1, Point p2);
            void drawSegment(Segment s);
            Prop getRegionFromId(uint32_t id);
            uint32_t getRegionId(Prop val);
        
            void plot();
    };

    // templated member functions must be defined in the header file
    template <typename Prop>
    Drawing<Prop>::Drawing() {

    }

    template <typename Prop>
    Drawing<Prop>::~Drawing() {
        
    }

    template <typename Prop>
    void Drawing<Prop>::drawRectangle(Point p1, Point p2) {
        drawSegment(p1, Point(p2.x, p1.y));
        drawSegment(Point(p2.x, p1.y), p2);
        drawSegment(p2, Point(p1.x, p2.y));
        drawSegment(Point(p1.x, p2.y), p1);
    }

    template <typename Prop>
    void Drawing<Prop>::drawRectangle(Point p1, double width, double height) {
        drawRectangle(p1, Point(p1.x + width, p1.y + height));
    }

    template <typename Prop>
    void Drawing<Prop>::drawCircle(Point p, double radius, uint32_t n_segments) {
        double angle = 0;
        double angle_step = 2 * PI / n_segments;
        Point p1 = Point(p.x + radius, p.y);
        Point p2;
        for (uint32_t i = 0; i < n_segments; i++) {
            angle += angle_step;
            p2 = Point(p.x + radius * cos(angle), p.y + radius * sin(angle));
            drawSegment(p1, p2);
            p1 = p2;
        }
    }

    template <typename Prop>
    void Drawing<Prop>::drawCircle(Circle c, uint32_t n_segments) {
        drawCircle(c.center, c.radius, n_segments);
    }

    template <typename Prop>
    void Drawing<Prop>::drawPolygon(Point* points, uint32_t n_points) {
        // check if the polygon self-intersects
        for (uint32_t i = 0; i < n_points; i++) {
            for (uint32_t j = i + 2; j < n_points; j++) {
                if (i == 0 && j == n_points - 1) {
                    continue;
                }
                if (Segment::segmentsIntersect(points[i], points[(i + 1) % n_points], points[j], points[(j + 1) % n_points])) {
                    nexit("Error: polygon self-intersects");
                }
            }
        }

        for (uint32_t i = 0; i < n_points - 1; i++) {
            drawSegment(points[i], points[i + 1]);
        }
        drawSegment(points[n_points - 1], points[0]);
    }

    template <typename Prop>
    void Drawing<Prop>::drawPolyLine(Point* points, uint32_t n_points) {
        for (uint32_t i = 0; i < n_points - 1; i++) {
            drawSegment(points[i], points[i + 1]);
        }
    }

    template <typename Prop>
    uint32_t Drawing<Prop>::getRegionId(Prop val) {
        if (region_map.find(val) == region_map.end()) {
            region_map[val] = region_map.size();
            // printf("region %f has id %lu\n", val, region_map[val]);
        }
        return region_map[val];
    }

    template <typename Prop>
    Prop Drawing<Prop>::getRegionFromId(uint32_t id) {
        for (auto it = region_map.begin(); it != region_map.end(); it++) {
            if (it->second == id) {
                return it->first;
            }
        }
        nexit("Error: region id not found");
        // unreachable
        exit(1);
    }

    template <typename Prop>
    void Drawing<Prop>::drawRegion(Point p, Prop val) {
        uint32_t region_id = getRegionId(val);
        regions.push_back(DrawingRegion(p, region_id));
    }

    template <typename Prop>
    void Drawing<Prop>::drawSegment(Point p1, Point p2) {
        // check if point is already in points
        bool p1_found = false;
        bool p2_found = false;
        uint32_t p1_id = 0;
        uint32_t p2_id = 0;
        for (uint32_t i = 0; i < points.size(); i++) {
            if (points[i] == p1) {
                p1_found = true;
                p1_id = i;
            }
            if (points[i] == p2) {
                p2_found = true;
                p2_id = i;
            }
        }
        if (!p1_found) {
            p1_id = points.size();
            points.push_back(p1);
        }
        if (!p2_found) {
            p2_id = points.size();
            points.push_back(p2);
        }
        
        segments.insert(DrawingSegment(p1_id, p2_id));
    }

    template <typename Prop>
    void Drawing<Prop>::drawSegment(Segment s) {
        drawSegment(s.p1, s.p2);
    }

    template <typename Prop>
    void Drawing<Prop>::plot() {
        // draw the mesh
        // returns zero on success else non-zero
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL");
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
        for (auto v : points) {
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
            SDL_SetRenderDrawColor(rend, 255, 0, 0, 255);
            for (DrawingSegment s : segments) {
                SDL_RenderDrawLine(rend, x_scale * points[s.p1].x + x_offset, y_scale * points[s.p1].y + y_offset, x_scale * points[s.p2].x + x_offset, y_scale * points[s.p2].y + y_offset);
            }
    

            // draw the points
            for (auto v : points) {
                SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
                // draw a square centered at the point
                SDL_Rect rect;
                rect.x = x_offset + v.x * x_scale - 2;
                rect.y = y_offset + v.y * y_scale - 2;
                rect.w = 4;
                rect.h = 4;
                SDL_RenderFillRect(rend, &rect);
            }

            // draw the regions
            for (DrawingRegion r : regions) {
                SDL_SetRenderDrawColor(rend, 0, 0, 255, 255);
                // draw a square centered at the point
                SDL_Rect rect;
                rect.x = x_offset + r.first.x * x_scale - 4;
                rect.y = y_offset + r.first.y * y_scale - 4;
                rect.w = 8;
                rect.h = 8;
                SDL_RenderFillRect(rend, &rect);
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

            // sleep for 1 second
            SDL_Delay(10);
        }
    }
}

#endif