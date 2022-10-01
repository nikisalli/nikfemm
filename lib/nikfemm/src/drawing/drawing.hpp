#ifndef NIK_DRAWING_HPP
#define NIK_DRAWING_HPP

#include <unordered_set>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#include "../geometry/point.hpp"
#include "../geometry/segment.hpp"
#include "drawing_region.hpp"
#include "drawing_segment.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"
#include "../mesh/vertex.hpp"

#define BOUNDARY_REGION PredefinedRegion(0)

namespace nikfemm {
    struct PredefinedRegion {
        double region_attribute;

        PredefinedRegion(double region_attribute) {
            this->region_attribute = region_attribute;
        }
    };

    struct Drawing {
        std::unordered_set<DrawingRegion> regions;
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
            void drawRegion(Point p, uint32_t region_attribute);
            void drawRegion(Point p, PredefinedRegion region);
            void drawSegment(Point p1, Point p2);
            void drawSegment(Segment s);
            void drawSegment(const Vertex &v1, const Vertex &v2);
            void drawSegment(const Vertex *v1, const Vertex *v2);
        
            void plot();
    };
}

#endif