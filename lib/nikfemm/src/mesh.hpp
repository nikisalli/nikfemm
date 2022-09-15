#ifndef NIK_MESH_HPP
#define NIK_MESH_HPP

#include <unordered_set>

#include "mesh_objects.hpp"


namespace nikfemm {
    struct Mesh {
        std::unordered_set<TriangleVertex*, std::hash<TriangleVertex*>, std::equal_to<TriangleVertex*>> vertices;
        std::unordered_set<TriangleElement*, std::hash<TriangleElement*>, std::equal_to<TriangleElement*>> elements;

        std::unordered_set<TriangleVertex*, std::hash<TriangleVertex*>, std::equal_to<TriangleVertex*>> boundary_vertices;

        Mesh();
        ~Mesh();

        void plot();
    };
}

#endif