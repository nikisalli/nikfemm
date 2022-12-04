#ifndef NIK_DRAWING_HPP
#define NIK_DRAWING_HPP

#include <unordered_set>
#include <utility>
#include <vector>
#include <unordered_map>
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
#include "../geometry/polygon.hpp"

namespace nikfemm {
    typedef std::pair<Point, uint32_t> DrawingRegion;

    template<typename Prop>
    struct Drawing {
        std::vector<Prop> region_map;
        std::vector<DrawingRegion> regions;
        std::vector<Point> points;
        std::vector<DrawingSegment> segments;
        std::vector<Polygon> polygons;

        Drawing();
        ~Drawing();

        /* drawing */
        public:
            void drawRectangle(Point p1, Point p2);
            void drawRectangle(Point p1, double width, double height);
            void drawCircle(Point p, double radius, uint32_t n_segments);
            void drawCircle(Circle c, uint32_t n_segments);
            void drawPolygon(Point* points, uint32_t n_points);
            void drawPolygon(const std::vector<Point>& points);
            void drawRegion(Point p, Prop val);
            const Prop* getRegionPtrFromId(uint32_t id) const;
            Prop getRegionFromId(uint32_t id) const;
            uint32_t getRegionId(Prop val);
            void addRefiningPoints();
            void plot();
        private:
            void drawSegment(Point p1, Point p2);
            void drawSegment(Segment s);
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
        polygons.push_back(Polygon({p1, Point(p2.x, p1.y), p2, Point(p1.x, p2.y)}));
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
        std::vector<Point> points;
        for (uint32_t i = 0; i < n_segments; i++) {
            angle += angle_step;
            p2 = Point(p.x + radius * cos(angle), p.y + radius * sin(angle));
            drawSegment(p1, p2);
            p1 = p2;
            points.push_back(p2);
        }
        polygons.push_back(Polygon(points));
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
        polygons.push_back(Polygon(points, n_points));
    }

    template <typename Prop>
    void Drawing<Prop>::drawPolygon(const std::vector<Point>& points) {
        drawPolygon((Point*) points.data(), points.size());
    }

    template <typename Prop>
    uint32_t Drawing<Prop>::getRegionId(Prop val) {
        for (uint32_t i = 0; i < region_map.size(); i++) {
            if (region_map[i] == val) {
                return i;
            }
        }
        region_map.push_back(val);
        return region_map.size() - 1;
    }

    template <typename Prop>
    Prop Drawing<Prop>::getRegionFromId(uint32_t id) const {
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
    const Prop* Drawing<Prop>::getRegionPtrFromId(uint32_t id) const {
        if (id >= region_map.size()) {
            nexit("Error: region id not found");
        }
        
        return &region_map.at(id);
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
        
        for (auto s : segments) {
            if (s.p1 == p1_id && s.p2 == p2_id) {
                printf("Warning: segment already exists\n");
                return;
            }
            if (Segment::segmentsIntersect(p1, p2, points[s.p1], points[s.p2]) && s.p1 != p1_id && s.p1 != p2_id && s.p2 != p1_id && s.p2 != p2_id) {
                nexit("Error: segment intersects another segment");
            }
        }
        segments.push_back(DrawingSegment(p1_id, p2_id));
    }

    template <typename Prop>
    void Drawing<Prop>::drawSegment(Segment s) {
        drawSegment(s.p1, s.p2);
    }

    template <typename Prop>
    void Drawing<Prop>::addRefiningPoints() {
        // find shortest segment first
        double length = std::numeric_limits<double>::max();
        for (auto it = segments.begin(); it != segments.end(); it++) {
            double len= Point::distance(points[it->p1], points[it->p2]);
            if (len < length) {
                length = len;
            }
        }

        length /= 100;

        uint32_t n_points = points.size();

        // for each corner, add as many points as possible and make sure they are at least 20 degrees apart
        // points should also stay away from the segments
        for (uint32_t i = 0; i < n_points; i++) {
            // find all segments that contain this point
            std::vector<DrawingSegment> segments_containing_point;
            for (auto it = segments.begin(); it != segments.end(); it++) {
                if (it->p1 == i || it->p2 == i) {
                    segments_containing_point.push_back(*it);
                }
            }

            // printf("point %u has %lu segments\n", i, segments_containing_point.size());

            for (uint32_t j = 0; j < 18; j++) {
                double angle = j * M_PI / 9;
                Point p = Point(points[i].x + length * cos(angle), points[i].y + length * sin(angle));
                
                // check distance to segments
                bool too_close = false;
                for (uint32_t k = 0; k < segments_containing_point.size(); k++) {
                    if (Segment::pointSegmentDistance(p, Segment(points[segments_containing_point[k].p1], points[segments_containing_point[k].p2])) < length / 2) {
                        too_close = true;
                        break;
                    }
                }
                if (!too_close) {
                    points.push_back(p);
                }
            }
        }
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
        float min_x = 1000000;
        float min_y = 1000000;
        float max_x = -1000000;
        float max_y = -1000000;
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
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * 2000 - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * 2000 - 0.5 * (max_y + min_y) * y_scale;

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