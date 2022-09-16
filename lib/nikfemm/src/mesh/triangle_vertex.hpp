#ifndef NIK_TRIANGLE_VERTEX_HPP
#define NIK_TRIANGLE_VERTEX_HPP

#include <cstdint>

#include "vertex.hpp"
#include "element.hpp"

namespace nikfemm {
    class TriangleVertex : public Vertex {
        public:

            Vertex* adjvert[18];
            Element* adjele[18];

            uint8_t adjvert_count = 0;
            uint8_t adjele_count = 0;

            TriangleVertex();
            TriangleVertex(Point p);
            TriangleVertex(double x, double y);
            ~TriangleVertex();

            void addAdjacentVertex(Vertex* v);
            void addAdjacentElement(Element* e);

            bool operator==(const Vertex& v) const;
            bool operator!=(const Vertex& v) const;
    };
}

namespace std {
    template <>
    struct hash<nikfemm::TriangleVertex> {
        inline std::size_t operator()(const nikfemm::TriangleVertex& v) const {
            return std::hash<nikfemm::Point>()(v.p);
        }
    };

    template <>
    struct hash<nikfemm::TriangleVertex*> {
        inline std::size_t operator()(const nikfemm::TriangleVertex* v) const {
            return std::hash<nikfemm::TriangleVertex>()(*v);
        }
    };

    template <>
    struct equal_to<nikfemm::TriangleVertex*> {
        inline bool operator()(const nikfemm::TriangleVertex* v1, const nikfemm::TriangleVertex* v2) const {
            if (v1 == v2) {
                return true;
            } else {
                return *v1 == *v2;
            }
        }
    };
}

#endif