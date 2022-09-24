#ifndef NIK_MESH_HPP
#define NIK_MESH_HPP

#include <unordered_set>
#include <vector>

#include "../drawing/drawing.hpp"
#include "vertex.hpp"
#include "../matrix/simple_vector.hpp"
#include "../matrix/csr.hpp"
#include "../matrix/coo.hpp"


namespace nikfemm {
    struct Mesh {
        Point center = Point(0, 0);
        double radius = 0;

        std::vector<Vertex*> vertices;
        std::vector<Vertex*> boundary_vertices;

        Mesh();
        ~Mesh();

        void plot();
        void mesh(Drawing &drawing);
        void addKelvinBoundaryConditions();
        void kelvinTransformCentered();
        void enumerateVertices();
        void getFemMatrix(MatCOO &coo);
        void getCoefficientVector(CV &b);
    };
}

#endif