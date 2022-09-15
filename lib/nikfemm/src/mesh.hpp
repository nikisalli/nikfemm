#ifndef NIK_MESH_HPP
#define NIK_MESH_HPP

#include <unordered_set>

#include "mesh_objects.hpp"


namespace nikfemm {
    struct Mesh {
        std::unordered_set<Vertex*, std::hash<Vertex*>, std::equal_to<Vertex*>> vertices;
        std::unordered_set<Element*, std::hash<Element*>, std::equal_to<Element*>> elements;

        std::unordered_set<Vertex*, std::hash<Vertex*>, std::equal_to<Vertex*>> boundary_vertices;
        std::unordered_set<Element*, std::hash<Element*>, std::equal_to<Element*>> boundary_elements;

        Mesh();
        ~Mesh();

        void plot();
    };
}

#endif