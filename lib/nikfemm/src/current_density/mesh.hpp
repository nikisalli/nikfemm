#ifndef NIK_CURRENT_DENSITYMESH_HPP
#define NIK_CURRENT_DENSITYMESH_HPP

#include <chrono>

#include "../triangle/triangle.h"

#include "../utils/utils.hpp"
#include "../drawing/drawing.hpp"
#include "../mesh/mesh.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/system.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct CurrentDensityMesh : Mesh<CurrentDensityProp> {
        CurrentDensityMesh();

        System<double> getFemSystem();  // tent function weights for energy minimization
        System<double> getFemSystemCotangentWeights();  // cotangent laplacian approximation weights
    };
}

#endif // NIK_CURRENT_DENSITYMESH_HPP