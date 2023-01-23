#ifndef NIK_CIRCLE_HPP
#define NIK_CIRCLE_HPP

#include <vector>
#include <cstdint>
#include <unordered_set>

#include "vector.hpp"
#include "common.hpp"

namespace nikfemm {
    struct Circle {
        public:
            Vector center;
            float radius;

            Circle(Vector center, float radius);
            Circle();

            bool operator==(const Circle& c) const;
            bool operator!=(const Circle& c) const;

            inline bool contains(Vector p) {
                return Vector::distance(center, p) <= radius;
            }
            static inline Circle getCircleFromPoints(Vector p1, Vector p2, Vector p3) {
                float bx = p2.x - p1.x;
                float by = p2.y - p1.y;
                float cx = p3.x - p1.x;
                float cy = p3.y - p1.y;

                float B = bx * bx + by * by;
                float C = cx * cx + cy * cy;
                float D = bx * cy - by * cx;
                Vector I = Vector((cy * B - by * C) / (2 * D), (bx * C - cx * B) / (2 * D));

                I.x += p1.x;
                I.y += p1.y;

                return Circle(I, Vector::distance(I, p1));
            }

            static inline Circle getCircleFromPoints(Vector p1, Vector p2) {
                Vector I = Vector((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
                return Circle(I, Vector::distance(I, p1));
            }

            inline bool containsPoints(std::vector<Vector> points) {
                for (auto p : points) {
                    if (!contains(p)) {
                        return false;
                    }
                }

                return true;
            }

            inline double circumference() {
                return 2 * PI * radius;
            }

            static Circle trivialCircleFromPoints(std::vector<Vector> points);

        protected:
            static Circle welzlHelper(std::vector<Vector> points, std::vector<Vector> R, uint32_t n);
        public:
            static Circle getMinimumEnclosingCircle(std::vector<Vector> points);
            static Circle getMinimumEnclosingCircle(std::unordered_set<Vector> points);
            static Circle getEnclosingCircle(std::vector<Vector> points);
    };
}

#endif