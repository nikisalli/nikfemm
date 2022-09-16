#include <assert.h>

#include "triangle_vertex.hpp"

namespace nikfemm {
    TriangleVertex::TriangleVertex() {

    }

    TriangleVertex::TriangleVertex(Point p) {
        p = p;
    }

    TriangleVertex::TriangleVertex(double x, double y) {
        p.x = x;
        p.y = y;
    }

    TriangleVertex::~TriangleVertex() {

    }

    void TriangleVertex::addAdjacentVertex(Vertex* v) {
        assert(adjvert_count < 18);
        // check if vertex is already in list
        for (int i = 0; i < adjvert_count; i++) {
            if (*adjvert[i] == *v) {
                return;
            }
        }
        adjvert[adjvert_count] = v;
        adjvert_count++;
    }

    void TriangleVertex::addAdjacentElement(Element* e) {
        assert(adjele_count < 18);
        // check if element is already in list
        for (int i = 0; i < adjele_count; i++) {
            if (*adjele[i] == *e) {
                return;
            }
        }
        adjele[adjele_count] = e;
        adjele_count++;
    }

    bool TriangleVertex::operator==(const Vertex& v) const {
        return p == v.p;
    }

    bool TriangleVertex::operator!=(const Vertex& v) const {
        return p != v.p;
    }
}