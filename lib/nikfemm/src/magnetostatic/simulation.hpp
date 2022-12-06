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
            std::vector<Vector> B;
            CV A;

            MagnetostaticSimulation();
            ~MagnetostaticSimulation();

            void solve();
        protected:
            void AplotRend(cv::Mat* image, double width, double height);
            void BplotRend(cv::Mat* image, double width, double height, bool plotMesh = false, bool plotRegions = false);
            void ScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
        public:
            void Aplot(uint32_t width, uint32_t height);
            void Bplot(uint32_t width, uint32_t height, bool plotMesh = false, bool plotRegions = false);
            void ScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void AplotToFile(uint32_t width, uint32_t height, std::string filename);
            void BplotToFile(uint32_t width, uint32_t height, std::string filename, bool plotMesh = false, bool plotRegions = false);
            void ScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh = false, bool plotRegions = false);
    };
}

#endif