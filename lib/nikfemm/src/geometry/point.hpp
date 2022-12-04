#ifndef NIK_POINT_HPP
#define NIK_POINT_HPP

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "../constants.hpp"

namespace nikfemm {
    struct Vector;
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
            Point operator*(const double d) const;
            Point operator/(const double d) const;
            
            // cast to vector
            operator Vector() const;

            static inline double double_oriented_area(const Point& p1, const Point& p2, const Point& p3) {
                return p1.x * p2.y - p1.x * p3.y - p2.x * p1.y + p2.x * p3.y + p3.x * p1.y - p3.x * p2.y;
            }

            static inline double distance(const Point p1, const Point p2) {
                return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
            }

            static inline double distance_squared(const Point p1, const Point p2) {
                return pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
            }
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