#ifndef NIK_MESH_OBJECTS_H
#define NIK_MESH_OBJECTS_H

#include <stdio.h>
#include <stdlib.h>
#include <cstdint>

#include <constants.hpp>

#include "geometry/point.hpp"
#include "mesh_objects.hpp"

#define EMPTY_ELEMENT nullptr
#define EMPTY_NODE nullptr
#define EMPTY_EDGE nullptr

namespace nikfemm {
    enum ElementType {
        ELEMENT_TYPE_TRIANGLE,
    };

    class Element;

    class Vertex {
        public:
            std::vector<Vertex*> adjvert;  // adjacent vertexes
            std::vector<Element*> adjele;  // adjacent elements

            /* properties */
            Point p;
            double A; // magnetic vector potential

            Vertex();
            Vertex(Point p);
            Vertex(double x, double y);
            ~Vertex();

            void addAdjacentVertex(Vertex* v);
            void addAdjacentElement(Element* e);

            static inline double distance(Vertex n1, Vertex n2) {
                return Point::distance(n1.p, n2.p);
            }

            bool operator==(const Vertex& v) const;
            bool operator!=(const Vertex& v) const;
    };

    class Element {
        public:
            double mu_r;

            ElementType type;

            std::vector<Element*> adjele;  // adjacent elements
            std::vector<Vertex*> vertices;   // vertices
        
            Element(Vertex* v1, Vertex* v2, Vertex* v3);
            Element(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4);
            Element(std::vector<Vertex*> vertices, std::vector<Element*> neighbors);
            Element();
            ~Element();

            void addAdjacentVertex(Vertex* v);
            void addAdjacentElement(Element* e);

            double getArea();

            bool operator==(const Element& e) const;
            bool operator!=(const Element& e) const;
    };

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
                return Vertex::distance(*e.v1, *e.v2);
            }

            bool operator==(const Edge& e) const;
            bool operator!=(const Edge& e) const;
            
    };
}

namespace std {
    /* hashers */
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

    template <>
    struct hash<nikfemm::Element> {
        inline std::size_t operator()(const nikfemm::Element& e) const {
            return std::hash<nikfemm::Vertex>()(*e.vertices[0]) ^ std::hash<nikfemm::Vertex>()(*e.vertices[1]) ^ std::hash<nikfemm::Vertex>()(*e.vertices[2]);
        }
    };

    template <>
    struct hash<nikfemm::Element*> {
        inline std::size_t operator()(const nikfemm::Element* e) const {
            return std::hash<nikfemm::Element>()(*e);
        }
    };

    /* comparators */

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

#endif // OBJECT_H