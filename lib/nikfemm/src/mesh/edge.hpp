#ifndef NIK_EDGE_HPP
#define NIK_EDGE_HPP

#define EMPTY_EDGE nullptr

#include <cstdint>

#include "vertex.hpp"
#include "element.hpp"

namespace nikfemm {
    class Edge {
        public:
            Vertex* v1;
            Vertex* v2;

            Element* neighbor1 = EMPTY_ELEMENT;
            Element* neighbor2 = EMPTY_ELEMENT;

            Edge(Vertex* v1, Vertex* v2);
            ~Edge();

            void addAdjacentElement(Element* e);
            uint32_t getAdjacentElementNum();
            bool hasVertices(Vertex* v1, Vertex* v2);

            static inline double length(Edge e) {
                return Point::distance(e.v1->p, e.v2->p);
            }

            bool operator==(const Edge& e) const;
            bool operator!=(const Edge& e) const;
            
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Edge> {
        inline std::size_t operator()(const nikfemm::Edge& e) const {
            return std::hash<nikfemm::Vertex>()(*e.v1) ^ std::hash<nikfemm::Vertex>()(*e.v2);
        }
    };

    template <>
    struct hash<nikfemm::Edge*> {
        inline std::size_t operator()(const nikfemm::Edge* e) const {
            return std::hash<nikfemm::Edge>()(*e);
        }
    };
}

#endif