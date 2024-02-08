#ifndef NIK_CURRENT_DENSITY_SIMULATION_HPP
#define NIK_CURRENT_DENSITY_SIMULATION_HPP

#include "../triangle/triangle.h"
#include "../../src/triangle/util.h"
#include "../constants.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"
#include "../drawing/drawing.hpp"
#include "../magnetostatic/mesh.hpp"
#include "mesh.hpp"

namespace nikfemm {
    class CurrentDensitySimulation {
        public:
            CurrentDensityMesh mesh;

            CurrentDensitySimulation(double depth);
            CurrentDensitySimulation();
            ~CurrentDensitySimulation();

            System<double> generateSystem(bool refine = true, double max_triangle_area = 1, int min_angle = 33);
            std::vector<double> solve(System<double>& system);
            void setVoltage(System<double>& system, Vector p, double V);
    };
}

#endif