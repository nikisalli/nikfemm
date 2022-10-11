#ifndef NIK_POINT_HPP
#define NIK_POINT_HPP

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "../constants.hpp"

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
            Point operator+(const Point& p) const;
            Point operator-(const Point& p) const;
            Point operator*(const double& d) const;
            Point operator/(const double& d) const;
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Point> {
        inline std::size_t operator()(const nikfemm::Point& p) const {
            return hash<double>()(p.x) ^ hash<double>()(p.y);
        }
    };

    template <>
    struct equal_to<nikfemm::Point> {
        inline bool operator()(const nikfemm::Point& p1, const nikfemm::Point& p2) const {
            return (p1.x - p2.x) < EPSILON && (p1.y - p2.y) < EPSILON;
        }
    };
}

#endif