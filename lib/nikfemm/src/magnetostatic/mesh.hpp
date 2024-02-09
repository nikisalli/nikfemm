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

#include "../triangle/triangle.h"

#include "../utils/utils.hpp"
#include "../drawing/drawing.hpp"
#include "../mesh/mesh.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/system.hpp"
#include "../algebra/non_linear.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticMesh : Mesh<MagnetostaticProp> {
        MagnetostaticMesh();

        System<NonLinearExpression> getFemSystem();
        void addDirichletInfiniteBoundaryConditions(System<NonLinearExpression>& system);
        void refineMeshAroundMagnets();
    };
}

#endif