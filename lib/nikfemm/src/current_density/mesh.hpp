#ifndef NIK_CURRENT_DENSITYMESH_HPP
#define NIK_CURRENT_DENSITYMESH_HPP

#include <chrono>
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
    struct CurrentDensityMesh : Mesh<CurrentDensityProp> {
        CurrentDensityMesh(double max_triangle_area);
        CurrentDensityMesh();
        ~CurrentDensityMesh();

        CurrentDensitySystem getFemSystem();  // tent function weights for energy minimization
        CurrentDensitySystem getFemSystemCotangentWeights();  // cotangent laplacian approximation weights
    };
}

#endif // NIK_CURRENT_DENSITYMESH_HPP