#ifndef NIK_DRAWING_HPP
#define NIK_DRAWING_HPP

#include <unordered_set>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#include "geometry/point.hpp"
#include "geometry/segment.hpp"
#include "drawing_region.hpp"
#include "drawing_segment.hpp"
#include "geometry/segment.hpp"
#include "geometry/circle.hpp"
#include "mesh/triangle_vertex.hpp"

#define BOUNDARY_REGION PredefinedRegion(-1)

namespace nikfemm {
    struct PredefinedRegion {
        int64_t region_id;

        PredefinedRegion(int64_t region_id) {
            this->region_id = region_id;
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
            void drawRegion(Point p, uint32_t region_id);
            void drawRegion(Point p, PredefinedRegion region);
            void drawSegment(Point p1, Point p2);
            void drawSegment(Segment s);
            void drawSegment(TriangleVertex v1, TriangleVertex v2);
            void drawSegment(TriangleVertex &v1, TriangleVertex &v2);
            void drawSegment(TriangleVertex *v1, TriangleVertex *v2);
        
            void plot();
    };
}

#endif