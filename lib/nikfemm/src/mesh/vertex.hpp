#ifndef NIK_VERTEX_HPP
#define NIK_VERTEX_HPP

#define EMPTY_NODE nullptr

#include <cstdint>
#include <algorithm>

#include "../geometry/point.hpp"
#include "../geometry/vector.hpp"

namespace nikfemm {
    template<typename Prop>
    class Vertex {
        public:
            /* properties */
            uint64_t id;
            Prop prop; // property of the region this vertex belongs to
            Point p; // position

            Vertex* adjvert[18];
            Prop adjprop[18];

            uint8_t adjvert_count = 0;
            uint8_t adjprop_count = 0;

            Vertex();
            Vertex(Point p);
            Vertex(double x, double y);
            ~Vertex();

            void addAdjacentVertex(Vertex* v);
            void addAdjacentProp(Prop prop);

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

    template<typename Prop>
    Vertex<Prop>::Vertex() {

    }

    template<typename Prop>
    Vertex<Prop>::Vertex(Point p) {
        p = p;
    }

    template<typename Prop>
    Vertex<Prop>::Vertex(double x, double y) {
        p.x = x;
        p.y = y;
    }

    template<typename Prop>
    Vertex<Prop>::~Vertex() {

    }

    template<typename Prop>
    void Vertex<Prop>::addAdjacentVertex(Vertex* v) {
        nassert(adjvert_count < 18, "Vertex::addAdjacentVertex: too many adjacent vertices");
        // check if vertex is already in list
        for (int i = 0; i < adjvert_count; i++) {
            if (*adjvert[i] == *v) {
                adjvert[i] = v;
                return;
            }
        }
        adjvert[adjvert_count] = v;
        adjvert_count++;
    }

    template<typename Prop>
    void Vertex<Prop>::addAdjacentProp(Prop prop) {
        nassert(adjprop_count < 18, "Vertex::addAdjacentProp: too many adjacent properties");
        // don't check if mu is already in list for performance
        adjprop[adjprop_count] = prop;
        adjprop_count++;
    }

    template<typename Prop>
    bool Vertex<Prop>::operator==(const Vertex& v) const {
        return p == v.p;
    }

    template<typename Prop>
    bool Vertex<Prop>::operator!=(const Vertex& v) const {
        return p != v.p;
    }

    template<typename Prop>
    double Vertex<Prop>::cellArea() {
        // sum of areas of triangles formed by vertex and adjacent vertices divided by 3
        double area = 0;
        // we know that the vertex is the center of the cell, use it to sort the adjacent vertices
        std::sort(adjvert, adjvert + adjvert_count, atanCompare(p));
        for (int i = 0; i < adjvert_count; i++) {
            area += geomArea(adjvert[i]->p, adjvert[(i + 1) % adjvert_count]->p, p);
        }
        return area / 3;
    }
}

#endif