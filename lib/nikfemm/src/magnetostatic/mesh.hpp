#ifndef NIK_MAGNETOSTATICMESH_HPP
#define NIK_MAGNETOSTATICMESH_HPP

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <math.h>
#include <set>
#include <opencv2/opencv.hpp>

extern "C" {
    #include "../../lib/triangle/triangle.h"
}

#include "../utils/utils.hpp"
#include "../drawing/drawing.hpp"
#include "../mesh/mesh.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/build_coo.hpp"
#include "algebra.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticMesh : Mesh<MagnetostaticProp> {
        MagnetostaticMesh(double max_triangle_area);
        MagnetostaticMesh();
        ~MagnetostaticMesh();

        MagnetostaticSystem getFemSystem();
        void addDirichletInfiniteBoundaryConditions(MagnetostaticSystem& system);
        void addDirichletZeroBoundaryConditions(MagnetostaticSystem& system, uint32_t id);
        void computeCurl(std::vector<Vector>& B, CV& A) const;
        void computeGrad(std::vector<Vector>& B, CV& A) const;
        void refineMeshAroundMagnets();
    };
}

#endif