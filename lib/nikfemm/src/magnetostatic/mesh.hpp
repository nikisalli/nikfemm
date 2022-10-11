#ifndef NIK_MESH_HPP
#define NIK_MESH_HPP

#include <unordered_set>
#include <vector>

#include "../drawing/drawing.hpp"
#include "vertex.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"
#include "../utils/utils.hpp"


namespace nikfemm {
    struct Mesh {
        Point center = Point(0, 0);
        double radius = 0;

        std::vector<Vertex*> vertices;
        std::vector<Vertex*> boundary_vertices;

        Mesh();
        ~Mesh();

        void plot();
        void Aplot();
        void Bplot();
        void mesh(Drawing &drawing);
        void addKelvinBoundaryConditions();
        void addDirichletBoundaryConditions(MatCOO &coo, CV &b);
        void kelvinTransformCentered();
        void enumerateVertices();
        void getFemSystem(MatCOO &coo, CV &b);
        void setField(CV &x);
        void computeCurl();
    };
}

#endif