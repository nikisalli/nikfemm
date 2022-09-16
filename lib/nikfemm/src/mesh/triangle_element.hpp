#ifndef NIK_TRIANGLE_ELEMENT_HPP
#define NIK_TRIANGLE_ELEMENT_HPP

#include "element.hpp"
#include "vertex.hpp"

namespace nikfemm {
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
}

#endif