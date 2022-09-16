#include "triangle_element.hpp"

namespace nikfemm {
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
}