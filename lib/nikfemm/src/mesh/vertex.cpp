#include <assert.h>
#include <algorithm>

#include "vertex.hpp"
#include "../geometry/geometry_common.hpp"

namespace nikfemm {
    Vertex::Vertex() {

    }

    Vertex::Vertex(Point p) {
        p = p;
    }

    Vertex::Vertex(double x, double y) {
        p.x = x;
        p.y = y;
    }

    Vertex::~Vertex() {

    }

    void Vertex::addAdjacentVertex(Vertex* v) {
        assert(adjvert_count < 18);
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

    void Vertex::addAdjacentMu(double mu) {
        assert(adjmuj_count < 18);
        // don't check if mu is already in list for performance
        adjmuj[adjmuj_count] = mu;
        adjmuj_count++;
    }

    bool Vertex::operator==(const Vertex& v) const {
        return p == v.p;
    }

    bool Vertex::operator!=(const Vertex& v) const {
        return p != v.p;
    }

    double Vertex::cellArea() {
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