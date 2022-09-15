#include <math.h>
#include <assert.h>

#include <mesh_objects.hpp>

namespace nikfemm {
    Vertex::~Vertex() {

    }

    void Vertex::addAdjacentVertex(Vertex* v) {

    }    

    void Vertex::addAdjacentElement(Element* e) {

    }

    bool Vertex::operator==(const Vertex& v) const {
        return false;
    }

    bool Vertex::operator!=(const Vertex& v) const {
        return false;
    }

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

    Element::~Element() {
        
    }

    double Element::getArea() {
        return 0;
    }

    bool Element::operator==(const Element& e) const {
        return false;
    }

    bool Element::operator!=(const Element& e) const {
        return true;
    }

    TriangleElement::TriangleElement(Vertex* v1, Vertex* v2, Vertex* v3) {
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;
    }
    
    TriangleElement::~TriangleElement() {
        
    }

    double TriangleElement::getArea() {
        double a = Point::distance(vertices[0]->p, vertices[1]->p);
        double b = Point::distance(vertices[1]->p, vertices[2]->p);
        double c = Point::distance(vertices[2]->p, vertices[0]->p);

        double s = (a + b + c) / 2;

        return sqrt(s * (s - a) * (s - b) * (s - c));
    }

    Point TriangleElement::getCenter() {
        double x = (vertices[0]->p.x + vertices[1]->p.x + vertices[2]->p.x) / 3;
        double y = (vertices[0]->p.y + vertices[1]->p.y + vertices[2]->p.y) / 3;

        return Point(x, y);
    }

    bool TriangleElement::operator==(const Element& e) const {
        if (e.type != ELEMENT_TYPE_TRIANGLE) {
            return false;
        }

        const TriangleElement& te = static_cast<const TriangleElement&>(e);

        // check same pointer
        if (this == &te) {
            return true;
        }

        // check if all vertices are the same
        for (Vertex* v : vertices) {
            bool found = false;
            for (Vertex* v2 : te.vertices) {
                if (v == v2) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }
        }

        return true;
    }

    bool TriangleElement::operator!=(const Element& e) const {
        return !(*this == e);
    }

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