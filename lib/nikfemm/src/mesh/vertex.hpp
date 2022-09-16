#ifndef NIK_VERTEX_HPP
#define NIK_VERTEX_HPP

#define EMPTY_NODE nullptr

#include "../geometry/point.hpp"
#include "element.hpp"

namespace nikfemm {
    class Vertex {
        public:
            /* properties */
            Point p;
            double A; // magnetic vector potential

            virtual ~Vertex() = 0;

            virtual void addAdjacentVertex(Vertex* v);
            virtual void addAdjacentElement(Element* e);

            virtual bool operator==(const Vertex& v) const;
            virtual bool operator!=(const Vertex& v) const;
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