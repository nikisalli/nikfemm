#ifndef NIK_MESH_HPP
#define NIK_MESH_HPP

#include <unordered_set>
#include <vector>

#include "../drawing/drawing.hpp"
#include "../mesh/vertex.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"
#include "../utils/utils.hpp"
#include "../constants.hpp"

namespace nikfemm {
    struct MagnetostaticProp {
        double mu; // permeability
        double J; // current density
        Vector M; // magnetization
        double A; // magnetic vector potential
        Vector B; // magnetic flux density

        // for sorting
        bool operator==(const MagnetostaticProp& other) const;
        bool operator!=(const MagnetostaticProp& other) const;
        bool operator<(const MagnetostaticProp& other) const;
    };

    // default material property
    const MagnetostaticProp vacuum_prop = {MU_0, 0, Vector(0, 0)};

    struct Mesh {
        Point center = Point(0, 0);
        double radius = 0;

        std::vector<Vertex<MagnetostaticProp>*> vertices;
        std::vector<Vertex<MagnetostaticProp>*> boundary_vertices;

        Mesh();
        ~Mesh();

        void plot();
        void Aplot();
        void Bplot();
        void mesh(Drawing<MagnetostaticProp> &drawing);
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