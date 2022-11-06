#ifndef NIK_CIRCLE_HPP
#define NIK_CIRCLE_HPP

#include <vector>
#include <cstdint>
#include <unordered_set>

#include "point.hpp"
#include "geometry_common.hpp"

namespace nikfemm {
    struct Circle {
        public:
            Point center;
            float radius;

            Circle(Point center, float radius);
            Circle();

            bool operator==(const Circle& c) const;
            bool operator!=(const Circle& c) const;

            inline bool contains(Point p) {
                return geomDistance(center, p) <= radius;
            }
            static inline Circle getCircleFromPoints(Point p1, Point p2, Point p3) {
                float bx = p2.x - p1.x;
                float by = p2.y - p1.y;
                float cx = p3.x - p1.x;
                float cy = p3.y - p1.y;

                float B = bx * bx + by * by;
                float C = cx * cx + cy * cy;
                float D = bx * cy - by * cx;
                Point I = Point((cy * B - by * C) / (2 * D), (bx * C - cx * B) / (2 * D));

                I.x += p1.x;
                I.y += p1.y;

                return Circle(I, geomDistance(I, p1));
            }

            static inline Circle getCircleFromPoints(Point p1, Point p2) {
                Point I = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
                return Circle(I, geomDistance(I, p1));
            }

            inline bool containsPoints(std::vector<Point> points) {
                for (auto p : points) {
                    if (!contains(p)) {
                        return false;
                    }
                }

                return true;
            }

            static Circle trivialCircleFromPoints(std::vector<Point> points);

        protected:
            static Circle welzlHelper(std::vector<Point> points, std::vector<Point> R, uint32_t n);
        public:
            static Circle getMinimumEnclosingCircle(std::vector<Point> points);
            static Circle getMinimumEnclosingCircle(std::unordered_set<Point> points);
    };
}

#endif