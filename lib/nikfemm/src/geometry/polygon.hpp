#ifndef NIK_POLYGON_HPP
#define NIK_POLYGON_HPP

#include <vector>

#include "point.hpp"

namespace nikfemm {
    class Polygon {
    public:
        std::vector<Point> points;

        Polygon();
        Polygon(const std::vector<Point>& points);
        Polygon(const Point* points, size_t n);
        Polygon(const Polygon& p);

        bool operator==(const Polygon& p) const;
        bool operator!=(const Polygon& p) const;

        bool contains(Point p) const;
        bool contains(Polygon p) const;
    };
}

#endif