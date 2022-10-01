#ifndef NIK_POINT_HPP
#define NIK_POINT_HPP

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <constants.hpp>

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
            Point operator+(const Point& p) const;
            Point operator-(const Point& p) const;
            Point operator*(const double& d) const;
            Point operator/(const double& d) const;

            static inline double distance(Point p1, Point p2) {
                return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
            }
            inline double distanceFromOrigin() {
                return sqrt(pow(x, 2) + pow(y, 2));
            }
            static inline double angle(Point a, Point b, Point c) {
                // safe angle calculation
                return fabs(atan2(c.y - b.y, c.x - b.x) - atan2(a.y - b.y, a.x - b.x));
            }
            static inline double area(Point a, Point b, Point c) {
                double ab = distance(a, b);
                double bc = distance(b, c);
                double ac = distance(a, c);
                double s = (ab + bc + ac) / 2;
                return sqrt(s * (s - ab) * (s - bc) * (s - ac));
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

    template <>
    struct equal_to<nikfemm::Point> {
        inline bool operator()(const nikfemm::Point& p1, const nikfemm::Point& p2) const {
            return (p1.x - p2.x) < EPSILON && (p1.y - p2.y) < EPSILON;
        }
    };
}

#endif