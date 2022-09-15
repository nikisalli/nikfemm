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
            /* properties */
            Point p;
            double A; // magnetic vector potential

            virtual ~Vertex() = 0;

            virtual void addAdjacentVertex(Vertex* v);
            virtual void addAdjacentElement(Element* e);

            virtual bool operator==(const Vertex& v) const;
            virtual bool operator!=(const Vertex& v) const;
    };

    class Element {
        public:
            double mu_r;

            ElementType type;

            virtual ~Element() = 0;

            virtual double getArea();

            virtual bool operator==(const Element& e) const;
            virtual bool operator!=(const Element& e) const;
    };

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

    class TriangleElement : public Element {
        public:
            double mu_r;

            ElementType type = ELEMENT_TYPE_TRIANGLE;

            Vertex* vertices[3];

            TriangleElement(Vertex* v1, Vertex* v2, Vertex* v3);
            ~TriangleElement();

            double getArea();
            Point getCenter();

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
                return Point::distance(e.v1->p, e.v2->p);
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

#endif // OBJECT_H