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
        CurrentDensityMesh(double max_triangle_area) {
            // default material property
            default_prop = {static_cast<float>(current_density_materials::copper)};
            this->max_triangle_area = max_triangle_area;
        }

        CurrentDensityMesh() {
            // default material property
            default_prop = {static_cast<float>(current_density_materials::copper)};
            max_triangle_area = 1e-0;
        }

        ~CurrentDensityMesh() {
            
        }

        CurrentDensitySystem getFemSystem();
        void addDirichletBoundaryConditions(CurrentDensitySystem& system, uint32_t id, double value);
    };

    CurrentDensitySystem CurrentDensityMesh::getFemSystem() {
        CurrentDensitySystem system = {
            BuildMatCOO<double>(data.numberofpoints),
            CV(data.numberofpoints)
        };

        auto start = std::chrono::high_resolution_clock::now();

        // COMPUTE FEM WEIGHTS
        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(data.numberofpoints);
        auto adjelems_props = std::vector<std::array<const CurrentDensityProp*, 18>>(data.numberofpoints);
        auto adjelems_count = std::vector<uint8_t>(data.numberofpoints, 0);

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = data.trianglelist[i][j];
                adjelems_ids[myid][adjelems_count[myid]] = i;
                adjelems_props[myid][adjelems_count[myid]++] = drawing.getRegionPtrFromId(data.triangleattributelist[i]);
            }
        }

        // surface integral
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            for (uint8_t j = 0; j < adjelems_count[i]; j++) {
                uint32_t v1, v2, v3;
                v1 = i;
                Elem myelem = data.trianglelist[adjelems_ids[i][j]];
                if (i == data.trianglelist[adjelems_ids[i][j]][0]) {
                    v2 = myelem[1];
                    v3 = myelem[2];
                } else if (i == myelem[1]) {
                    v2 = myelem[2];
                    v3 = myelem[0];
                } else if (i == myelem[2]) {
                    v2 = myelem[0];
                    v3 = myelem[1];
                } else {
                    nexit("error: vertex not found in element");
                }

                double oriented_area = data.getDoubleOrientedArea(v1, v2, v3);

                if (oriented_area < 0) {
                    std::swap(v2, v3);
                }

                double double_area = data.getDoubleOrientedArea(v1, v2, v3);

                // elements are only added if they are in the upper triangle because the matrix is symmetric and this saves half the memory
                double b1 = (data.pointlist[v2].y - data.pointlist[v3].y) / double_area;
                double c1 = (data.pointlist[v3].x - data.pointlist[v2].x) / double_area;
                double b2 = (data.pointlist[v3].y - data.pointlist[v1].y) / double_area;
                double c2 = (data.pointlist[v1].x - data.pointlist[v3].x) / double_area;
                double b3 = (data.pointlist[v1].y - data.pointlist[v2].y) / double_area;
                double c3 = (data.pointlist[v2].x - data.pointlist[v1].x) / double_area;

                if (v1 >= i) system.A(i, v1) += double_area * (b1 * b1 + c1 * c1) * 0.5 * (-adjelems_props[i][j]->sigma);
                if (v2 >= i) system.A(i, v2) += double_area * (b2 * b1 + c2 * c1) * 0.5 * (-adjelems_props[i][j]->sigma);
                if (v3 >= i) system.A(i, v3) += double_area * (b3 * b1 + c3 * c1) * 0.5 * (-adjelems_props[i][j]->sigma);

                // set the b vector
            }
            system.b.val[i] = 0;
        }

        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("FEM matrix construction took %d ms", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        return system;
    }
}

#endif // NIK_CURRENT_DENSITYMESH_HPP