#ifndef NIK_SIMULATION_H
#define NIK_SIMULATION_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "../lib/triangle/triangle.h"
#include "../src/triangle/util.h"

#include <constants.hpp>

#include "../src/geometry/segment.hpp"
#include "../src/geometry/circle.hpp"
#include "../src/drawing/drawing.hpp"
#include "../src/mesh/mesh.hpp"

namespace nikfemm {
    class Simulation {
        protected:
            Mesh mesh;
        
        public:
            Simulation();
            ~Simulation();

            /* meshing */
            void generateMesh(Drawing drawing);
    };
}

#endif