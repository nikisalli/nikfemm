#ifndef NIK_CURRENT_DENSITYMESH_HPP
#define NIK_CURRENT_DENSITYMESH_HPP

#include <chrono>

#ifdef NIKFEMM_USE_OPENCV
#include <opencv2/opencv.hpp>
#endif

#include "../triangle/triangle.h"

#include "../utils/utils.hpp"
#include "../drawing/drawing.hpp"
#include "../mesh/mesh.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/system.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct CurrentDensityMesh : Mesh<CurrentDensityProp> {
        double depth = 1;
        CurrentDensityProp default_prop = {static_cast<float>(current_density_materials::copper)};

        System<double> getFemSystem();  // tent function weights for energy minimization
        System<double> getFemSystemCotangentWeights();  // cotangent laplacian approximation weights
    };
}

#endif // NIK_CURRENT_DENSITYMESH_HPP