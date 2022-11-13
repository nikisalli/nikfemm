#ifndef NIK_SIMULATION_H
#define NIK_SIMULATION_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "../../lib/triangle/triangle.h"
#include "../../src/triangle/util.h"
#include "../constants.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"
#include "../drawing/drawing.hpp"
#include "../magnetostatic/mesh.hpp"
#include "mesh.hpp"

namespace nikfemm {
    class MagnetostaticSimulation {
        protected:
        
        public:
            MagnetostaticMesh mesh;
            MagnetostaticSimulation();
            ~MagnetostaticSimulation();

            /* meshing */
            void generateMesh();
    };
}

#endif