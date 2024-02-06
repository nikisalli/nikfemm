#ifndef NIK_MAGNETOSTATIC_SIMULATION_H
#define NIK_MAGNETOSTATIC_SIMULATION_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "../triangle/triangle.h"

#include "../../src/triangle/util.h"
#include "../constants.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"
#include "../drawing/drawing.hpp"
#include "../magnetostatic/mesh.hpp"
#include "mesh.hpp"

namespace nikfemm {
    struct SurroundingRegionBlockIntegralAssets {
        std::vector<std::array<uint32_t, 18>> adjelems_ids;
        std::vector<uint8_t> adjelems_count;
        std::vector<double> elem_field_errors;
        double min_err;
        double max_err;
        std::vector<bool> vertex_inside_integration_region;
        std::vector<bool> vertex_inside_integration_region_with_boundary;
        std::vector<bool> vertex_inside_boundary_region;
        std::vector<bool> vertex_inside_boundary_region_hole;
    };

    struct StressTensor {
        Vector Force;
        double Torque;
    };

    class MagnetostaticSimulation {
        public:
            MagnetostaticMesh mesh;
            std::vector<Vector> B;
            CV A;
            double depth = 1.0;

            MagnetostaticSystem generateSystem(bool refine = true, double max_triangle_area = 1, int min_angle = 33);
            void solve(MagnetostaticSystem& system);
            Vector computeForceIntegrals(Vector p);
            double computeTorqueIntegral(Vector p, Vector center);
            StressTensor computeStressIntegral(Vector p, Vector center);
            Vector computeBarycenter(Vector p);
            double computeInertiaMoment(Vector p, Vector center, double density);
            double computeMass(Vector p, double density);
        protected:
            double computeAreaIntegral(Vector p);
            Vector computeForceIntegrals(SurroundingRegionBlockIntegralAssets assets, Vector p);
            double computeTorqueIntegral(SurroundingRegionBlockIntegralAssets assets, Vector p, Vector center);
            Polygon getInnermostPolygon(Vector p);
#ifdef NIKFEMM_USE_OPENCV
            void AplotRend(cv::Mat* image, double width, double height);
            void BplotRend(cv::Mat* image, double width, double height, bool plotMesh = false, bool plotRegions = false, double maxB = NAN, double minB = NAN, bool plotCurves = false, uint32_t curve_number = 100);
            void ElemScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void NodeScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
#endif
            static void updateMu(std::vector<const MagnetostaticProp*>& props, std::vector<float>& mu, std::vector<Vector>& B, double residual, uint32_t iter);
            auto getSurroundingRegionBlockIntegralAssets(Vector p);
        public:
#ifdef NIKFEMM_USE_OPENCV
            void Aplot(uint32_t width, uint32_t height);
            void Bplot(uint32_t width, uint32_t height, bool plotMesh = false, bool plotRegions = false, double maxB = NAN, double minB = NAN, bool waitkey = true, bool plotCurves = false, uint32_t curve_number = 100);
            void ElemScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void NodeScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh = false, bool plotRegions = false);
            void AplotToFile(uint32_t width, uint32_t height, std::string filename);
            void BplotToFile(uint32_t width, uint32_t height, std::string filename, bool plotMesh = false, bool plotRegions = false, double maxB = NAN, double minB = NAN, bool plotCurves = false, uint32_t curve_number = 100);
            void ElemScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh = false, bool plotRegions = false);
            void NodeScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh = false, bool plotRegions = false);
#endif
    };
}

#endif