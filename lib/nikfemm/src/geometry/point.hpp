#ifndef NIK_POINT_HPP
#define NIK_POINT_HPP

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "geometry_common.hpp"

namespace nikfemm {
    class Point {
        public:
            double x;
            double y;

            Point(double x, double y);
            Point();
            ~Point();

            bool operator==(const Point& p) const;
            bool operator!=(const Point& p) const;
            bool operator<(const Point& p) const;

            static double distance(Point p1, Point p2) {
                return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
            }
            double distanceFromOrigin() {
                return sqrt(pow(x, 2) + pow(y, 2));
            }
            static Orientation orientation(Point p1, Point p2, Point p3);
            std::vector<Point> getWelzlPoints(std::vector<Point> points);
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Point> {
        inline std::size_t operator()(const nikfemm::Point& p) const {
            return hash<double>()(p.x) ^ hash<double>()(p.y);
        }
    };
}

#endif