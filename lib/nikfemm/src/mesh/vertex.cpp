#include <assert.h>

#include "vertex.hpp"

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
                return;
            }
        }
        adjvert[adjvert_count] = v;
        adjvert_count++;
    }

    void Vertex::addAdjacentMu(double mu) {
        assert(adjmu_r_count < 18);
        // don't check if mu is already in list for performance
        adjmu_r[adjmu_r_count] = mu;
        adjmu_r_count++;
    }

    bool Vertex::operator==(const Vertex& v) const {
        return p == v.p;
    }

    bool Vertex::operator!=(const Vertex& v) const {
        return p != v.p;
    }
}