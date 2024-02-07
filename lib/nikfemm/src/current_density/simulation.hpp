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
            std::vector<double> V;

            CurrentDensitySimulation(double depth);
            CurrentDensitySimulation();
            ~CurrentDensitySimulation();

            System<double> generateSystem(bool refine = true, double max_triangle_area = 1, int min_angle = 33);
            void solve(System<double>& system);
            void setVoltage(System<double>& system, Vector p, double V);
        protected:
#ifdef NIKFEMM_USE_OPENCV
            void VplotRend(cv::Mat* image, double width, double height);
            void JplotRend(cv::Mat* image, double width, double height);
            void EplotRend(cv::Mat* image, double width, double height, bool plotMesh = false, bool plotRegions = false, double maxB = NAN, double minB = NAN, bool plotCurves = false, uint32_t curve_number = 100);
            void ElemScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void NodeScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
#endif
        public:
#ifdef NIKFEMM_USE_OPENCV
            void Vplot(uint32_t width, uint32_t height);
            void Jplot(uint32_t width, uint32_t height);
            void Eplot(uint32_t width, uint32_t height, bool plotMesh = false, bool plotRegions = false, double maxB = NAN, double minB = NAN, bool waitkey = true, bool plotCurves = false, uint32_t curve_number = 100);
            void ElemScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void NodeScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void VplotToFile(uint32_t width, uint32_t height, std::string filename);
            void JplotToFile(uint32_t width, uint32_t height, std::string filename);
            void EplotToFile(uint32_t width, uint32_t height, std::string filename, bool plotMesh = false, bool plotRegions = false, double maxB = NAN, double minB = NAN, bool plotCurves = false, uint32_t curve_number = 100);
            void ElemScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh = false, bool plotRegions = false);
            void NodeScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh = false, bool plotRegions = false);
#endif
    };
}

#endif