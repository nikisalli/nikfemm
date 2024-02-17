#ifndef NIK_CURRENT_DENSITY_MULTILAYER_SIMULATION_HPP
#define NIK_CURRENT_DENSITY_MULTILAYER_SIMULATION_HPP

#include <vector>

#include "mesh.hpp"

namespace nikfemm {
    struct CurrentDensityInterconnection {
        // a connection between two points in two layers with a resistance
        Vector p1;
        Vector p2;
        uint64_t layer1_id; // corresponds to the index of the mesh in the meshes vector or system in the systems vector or depth in the depths vector
        uint64_t layer2_id; // same as above
        double R; // resistance
    };

    class MultiLayerCurrentDensitySimulation {
        public:
            std::vector<CurrentDensityMesh> meshes;
            std::vector<CurrentDensityInterconnection> interconnections;

            MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths);
            MultiLayerCurrentDensitySimulation(uint32_t num_layers);

            System<double> generateSystem(bool refine = true, double max_triangle_area = 1, int min_angle = 33);
            std::vector<double> solve(System<double>& system);
            void setVoltage(System<double>& system, Vector p, double V, uint64_t layer_id);
            std::vector<double> computePowerDensity(std::vector<double>& V);
    };
}

#endif