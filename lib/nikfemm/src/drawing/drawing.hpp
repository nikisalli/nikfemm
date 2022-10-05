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

#include "../geometry/point.hpp"
#include "../geometry/segment.hpp"
#include "drawing_segment.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"
#include "../mesh/vertex.hpp"

#define BOUNDARY_REGION 0.0f

namespace nikfemm {
    typedef std::pair<Point, uint64_t> DrawingRegion;

    struct Drawing {
        std::map<double, uint64_t> region_map;
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
            void drawRegion(Point p, double val);
            void drawSegment(Point p1, Point p2);
            void drawSegment(Segment s);
            void drawSegment(const Vertex &v1, const Vertex &v2);
            void drawSegment(const Vertex *v1, const Vertex *v2);
            double getRegionFromId(uint64_t id);
            uint64_t getRegionId(double val);
        
            void plot();
    };
}

#endif