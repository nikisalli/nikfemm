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
        double R;
    };

    class MultiLayerCurrentDensitySimulation {
        public:
            std::vector<CurrentDensityMesh> meshes;
            CV V;
            std::vector<double> depths;
            std::vector<CurrentDensityInterconnection> interconnections;

            MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths, std::vector<double> max_triangle_areas);
            MultiLayerCurrentDensitySimulation(uint32_t num_layers);
            ~MultiLayerCurrentDensitySimulation();

            CurrentDensitySystem generateSystem(bool refine = true);
            void solve(CurrentDensitySystem& system);
            void setVoltage(CurrentDensitySystem& system, Vector p, double V, uint64_t layer_id);
        protected:
            void VplotRend(cv::Mat* image, double width, double height, uint64_t layer_id, double maxV, double minV);
        public:
            void Vplot(uint32_t width, uint32_t height);
    };
}

#endif