#include <math.h>
#include <assert.h>

#include <mesh_objects.hpp>

namespace nikfemm {
    Vertex::Vertex() {
        p = Point(0, 0);
    }

    Vertex::Vertex(Point p) {
        this->p = p;
    }

    Vertex::Vertex(double x, double y) {
        p = Point(x, y);
    }

    Vertex::~Vertex() {
        
    }

    void Vertex::addAdjacentVertex(Vertex* v) {
        // check if v is already in adjvert
        for (Vertex* vert : adjvert) {
            if (vert == v) {
                return;
            }
        }
        adjvert.push_back(v);
    }

    void Vertex::addAdjacentElement(Element* e) {
        // check if e is already in adjele
        for (Element* ele : adjele) {
            if (ele == e) {
                return;
            }
        }
        adjele.push_back(e);
    }

    bool Vertex::operator==(const Vertex& v) const {
        return p == v.p;
    }

    bool Vertex::operator!=(const Vertex& v) const {
        return p != v.p;
    }

    Element::Element(Vertex* v1, Vertex* v2, Vertex* v3) {
        this->vertices.push_back(v1);
        this->vertices.push_back(v2);
        this->vertices.push_back(v3);
    }

    Element::Element(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) {
        this->vertices.push_back(v1);
        this->vertices.push_back(v2);
        this->vertices.push_back(v3);
        this->vertices.push_back(v4);
    }

    Element::Element(std::vector<Vertex*> vertices, std::vector<Element*> neighbors) {
        this->adjele = neighbors;
        this->vertices = vertices;
    }

    Element::Element() {
        adjele = std::vector<Element*>();
        vertices = std::vector<Vertex*>();
    }

    Element::~Element() {
        
    }

    double Element::getArea() {
        double a = Vertex::distance(*vertices[0], *vertices[1]);
        double b = Vertex::distance(*vertices[1], *vertices[2]);
        double c = Vertex::distance(*vertices[2], *vertices[0]);

        double s = (a + b + c) / 2.0;
        return sqrt(s * (s - a) * (s - b) * (s - c));
    }

    void Element::addAdjacentVertex(Vertex* v) {
        // check if v is already in vertices
        for (Vertex* vert : vertices) {
            if (vert == v) {
                return;
            }
        }
        vertices.push_back(v);
    }

    void Element::addAdjacentElement(Element* e) {
        // check if e is already in adj
        for (Element* ele : adjele) {
            if (ele == e) {
                return;
            }
        }
        adjele.push_back(e);
    }

    bool Element::operator==(const Element& e) const {
        return vertices == e.vertices;
    }

    bool Element::operator!=(const Element& e) const {
        return vertices != e.vertices;
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