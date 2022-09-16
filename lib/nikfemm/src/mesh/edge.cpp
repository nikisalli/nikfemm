#include <assert.h>

#include "edge.hpp"

namespace nikfemm {
    Edge::Edge(Vertex* v1, Vertex* v2) {
        this->v1 = v1;
        this->v2 = v2;
        this->neighbor1 = EMPTY_ELEMENT;
        this->neighbor2 = EMPTY_ELEMENT;
    }

    Edge::~Edge() {
        
    }

    void Edge::addAdjacentElement(Element* e) {
        assert(!(neighbor1 != EMPTY_ELEMENT && neighbor2 != EMPTY_ELEMENT));
        // check if e is already in neighbors
        if (neighbor1 == e || neighbor2 == e) {
            return;
        }
        if (neighbor1 == EMPTY_ELEMENT) {
            neighbor1 = e;
        } else if (neighbor2 == EMPTY_ELEMENT) {
            neighbor2 = e;
        }
    }

    uint32_t Edge::getAdjacentElementNum() {
        uint32_t num = 0;
        if (neighbor1 != EMPTY_ELEMENT) {
            num++;
        }
        if (neighbor2 != EMPTY_ELEMENT) {
            num++;
        }
        return num;
    }

    bool Edge::operator==(const Edge& e) const {
        return (*v1 == *(e.v1) && *v2 == *(e.v2)) || (*v1 == *(e.v2) && *v2 == *(e.v1));
    }

    bool Edge::operator!=(const Edge& e) const {
        return !(*this == e);
    }

    bool Edge::hasVertices(Vertex* v1, Vertex* v2) {
        return (*this->v1 == *v1 && *this->v2 == *v2) || (*this->v1 == *v2 && *this->v2 == *v1);
    }
}