#ifndef NIK_VERTEX_HPP
#define NIK_VERTEX_HPP

#define EMPTY_NODE nullptr

#include <cstdint>

#include "../geometry/point.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    class Vertex {
        public:
            /* properties */
            uint64_t id;
            double muj; // relative permeability * current density
            double A; // magnetic vector potential
            vec2 B; // magnetic flux density
            Point p;

            Vertex* adjvert[18];
            double adjmuj[18];

            uint8_t adjvert_count = 0;
            uint8_t adjmuj_count = 0;

            Vertex();
            Vertex(Point p);
            Vertex(double x, double y);
            ~Vertex();

            void addAdjacentVertex(Vertex* v);
            void addAdjacentMu(double mu);

            bool operator==(const Vertex& v) const;
            bool operator!=(const Vertex& v) const;

            double cellArea();

            // comparison operator for sorting based on atan2
            struct atanCompare {
                atanCompare(Point center) {
                    this->center = center;
                }

                Point center;
                bool operator()(const Vertex* v1, const Vertex* v2) const {
                    return atan2(v1->p.y - center.y, v1->p.x - center.x) < atan2(v2->p.y - center.y, v2->p.x - center.x);
                }
            };
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Vertex> {
        inline std::size_t operator()(const nikfemm::Vertex& v) const {
            return std::hash<nikfemm::Point>()(v.p);
        }
    };

    template <>
    struct hash<nikfemm::Vertex*> {
        inline std::size_t operator()(const nikfemm::Vertex* v) const {
            return std::hash<nikfemm::Vertex>()(*v);
        }
    };

    template <>
    struct equal_to<nikfemm::Vertex*> {
        inline bool operator()(const nikfemm::Vertex* v1, const nikfemm::Vertex* v2) const {
            if (v1 == v2) {
                return true;
            } else {
                return *v1 == *v2;
            }
        }
    };
}

#endif