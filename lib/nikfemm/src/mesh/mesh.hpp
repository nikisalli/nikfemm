#ifndef NIK_MESH_HPP
#define NIK_MESH_HPP

#include <unordered_set>
#include <vector>

#include "../drawing/drawing.hpp"
#include "triangle_element.hpp"
#include "triangle_vertex.hpp"


namespace nikfemm {
    struct Mesh {
        Point center = Point(0, 0);
        double radius = 0;

        std::unordered_set<TriangleVertex*, std::hash<TriangleVertex*>, std::equal_to<TriangleVertex*>> vertices;
        std::unordered_set<TriangleElement*, std::hash<TriangleElement*>, std::equal_to<TriangleElement*>> elements;

        std::vector<TriangleVertex*> boundary_vertices;

        Mesh();
        ~Mesh();

        void plot();
        void mesh(Drawing &drawing);
        void addKelvinBoundaryConditions();
        void kelvinTransformCentered();
    };
}

#endif