#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <set>

#include <omp.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "../constants.hpp"
#include "multilayer_simulation.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/vector.hpp"
#include "../algebra/build_coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths, std::vector<double> max_triangle_areas) {
        if (depths.size() != num_layers) {
            throw std::invalid_argument("depths.size() != num_layers");
        }
        if (max_triangle_areas.size() != num_layers) {
            throw std::invalid_argument("max_triangle_areas.size() != num_layers");
        }

        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].max_triangle_area = max_triangle_areas[i];
            this->depths.push_back(depths[i]);
        }
    }

    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers) {
        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].max_triangle_area = 1;
            this->depths.push_back(0);
        }
    }

    MultiLayerCurrentDensitySimulation::~MultiLayerCurrentDensitySimulation() {

    }

    CurrentDensitySystem MultiLayerCurrentDensitySimulation::generateSystem(bool refine) {
        if (refine) {
            for (auto& mesh : meshes) {
                mesh.drawing.addRefiningPoints();
            }
        }

        // add points where the interconnections are so that we can have the precise point in the mesh for it
        // insteaad of relying on the nearest point
        for (auto& interconnection : interconnections) {
            // the drawing algorithm will check internally if the point is already in the mesh
            // I could check if the point is inside the simulation domain but it will error out anyways if it is not
            // so I will just let it be
            meshes[interconnection.layer1_id].drawing.drawPoint(interconnection.p1);
            meshes[interconnection.layer2_id].drawing.drawPoint(interconnection.p2);
        }

        for (uint64_t i = 0; i < meshes.size(); i++) {
            meshes[i].mesh();
            meshes[i].computeEpsilon();
            nloginfo("Meshed layer %d with %d nodes and %d elements", i, meshes[i].data.numberofpoints, meshes[i].data.numberoftriangles);
        }

        std::vector<CurrentDensitySystem> systems;

        for (uint64_t i = 0; i < meshes.size(); i++) {
            systems.push_back(meshes[i].getFemSystem());
        }

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& system : systems) {
            offsets.push_back(offset);
            offset += system.b.val.size();
        }

        printf("offsets: ");
        for (auto& offset : offsets) {
            printf("%d ", offset);
        }
        printf("\n");

        fflush(stdout);

        // now we need to merge the systems together using the interconnections
        struct interconnection_indices {
            uint64_t mesh1_id;
            uint64_t mesh2_id;
            uint64_t node1_id;
            uint64_t node2_id;
            double R;
        };
        std::vector<interconnection_indices> interconnection_indices_list;

        for (auto& interconnection : interconnections) {
            // find the nearest node in the mesh
            uint64_t node1_id = 0;
            double node1_dist = INFINITY;
            uint64_t node2_id = 0;
            double node2_dist = INFINITY;

            for (uint64_t i = 0; i < meshes[interconnection.layer1_id].data.numberofpoints; i++) {
                double dist = 0;

                dist = (meshes[interconnection.layer1_id].data.pointlist[i] - interconnection.p1).norm();
                if (dist < node1_dist) {
                    node1_id = i;
                    node1_dist = dist;
                }

                dist = (meshes[interconnection.layer1_id].data.pointlist[i] - interconnection.p2).norm();
                if (dist < node2_dist) {
                    node2_id = i;
                    node2_dist = dist;
                }
            }

            nloginfo("Nearest node to p1 is %d", node1_id);
            nloginfo("Nearest node to p2 is %d", node2_id);

            interconnection_indices_list.push_back({interconnection.layer1_id, interconnection.layer2_id, node1_id, node2_id, interconnection.R});
        }

        // preallocate the size of the merged system by summing the sizes of the systems
        uint64_t merged_system_nnz = 0;
        uint64_t merged_system_size = 0;
        printf("merged_system_nnz is: %d\n", merged_system_nnz);
        printf("merged_system_size is: %d\n", merged_system_size);
        for (auto& system : systems) {
            merged_system_nnz += system.A.elems.size();
            merged_system_size += system.b.val.size();
        }

        CurrentDensitySystem merged_system = {
            BuildMatCOO<double>(merged_system_size),
            CV(merged_system_size)
        };

        printf("merged_system_nnz: %d\n", merged_system_nnz);
        printf("merged_system_size: %d\n", merged_system_size);
        merged_system.A.elems.reserve(merged_system_nnz);
        merged_system.b.val.reserve(merged_system_size);

        // merge the systems together like a block matrix
        offset = 0; // reusing offset variable

        for (uint64_t i = 0; i < meshes.size(); i++) {
            CurrentDensitySystem& system = systems[i];
            CurrentDensityMesh& mesh = meshes[i];

            uint64_t system_size = mesh.data.numberofpoints;

            // add the system to the merged system
            // A
            for (auto& elem : system.A.elems) {
                uint64_t m = elem.first >> 32;
                uint64_t n = elem.first & 0xFFFFFFFF;

                // printf("adding element (%d, %d) of layer %d of value %.17g to (%d, %d) of merged system\n", m, n, i, elem.second, m + offset, n + offset);

                merged_system.A(m + offset, n + offset) = elem.second;
            }

            // b
            for (uint64_t j = 0; j < system_size; j++) {
                merged_system.b.val.push_back(system.b.val[j]);
            }

            offset += system_size;
        }

        printf("b.size(): %d\n", merged_system.b.val.size());

        // add the interconnections to the merged system
        
        for (auto& interconnection : interconnection_indices_list) {
            // add the interconnection to the merged system
            uint64_t m1 = interconnection.node1_id + offsets[interconnection.mesh1_id];
            uint64_t m2 = interconnection.node2_id + offsets[interconnection.mesh2_id];
            double R = interconnection.R;

            // get all the elements with m as the row
            for (auto& elem : merged_system.A.elems) {
                uint64_t elem_m = elem.first >> 32;
                uint64_t elem_n = elem.first & 0xFFFFFFFF;

                if (elem_m == m1) {
                    if (elem_n == m1) {
                        printf("(%d, %d) of layer %d of value %.17g is being changed to %.17g\n", elem_m, elem_n, interconnection.mesh1_id, elem.second, (R * elem.second) - 1);
                        elem.second = (R * elem.second) +1;
                    } else {
                        printf("(%d, %d) of layer %d of value %.17g is being changed to %.17g\n", elem_m, elem_n, interconnection.mesh1_id, elem.second, R * elem.second);
                        elem.second *= R;
                    }
                } else if (elem_m == m2) {
                    if (elem_n == m2) {
                        printf("(%d, %d) of layer %d of value %.17g is being changed to %.17g\n", elem_m, elem_n, interconnection.mesh2_id, elem.second, (R * elem.second) - 1);
                        elem.second = (R * elem.second) +1;
                    } else {
                        printf("(%d, %d) of layer %d of value %.17g is being changed to %.17g\n", elem_m, elem_n, interconnection.mesh2_id, elem.second, R * elem.second);
                        elem.second *= R;
                    }
                }                        
            }

            // now subtract the corresponding element from the other system
            if (m2 >= m1) {
                merged_system.A(m1, m2) = -1;
                // printf("adding separate element (%d, %d) of layer %d of value %.17g to (%d, %d) of merged system\n", m1, m2, interconnection.mesh1_id, -1, m1, m2);
            } else {
                merged_system.A(m2, m1) = -1;
                // printf("adding separate element (%d, %d) of layer %d of value %.17g to (%d, %d) of merged system\n", m2, m1, interconnection.mesh2_id, -1, m2, m1);
            }
        }

        for (auto& elem : merged_system.A.elems) {
            printf("merged_system.A(%d, %d) = %.17g\n", elem.first >> 32, elem.first & 0xFFFFFFFF, elem.second);
        }       

        nloginfo("generated merged system with %d unknowns", merged_system.b.val.size());

        return merged_system;
    }

    void MultiLayerCurrentDensitySimulation::setVoltage(CurrentDensitySystem& system, Vector p, double V, uint64_t layer_id) {
        // find the closest node
        int32_t closest_node = -1;
        double closest_distance = INFINITY;
        auto& mesh = meshes[layer_id];

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& mesh : meshes) {
            offsets.push_back(offset);
            offset += mesh.data.numberofpoints;
        }

        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            double distance = (mesh.data.pointlist[i].x - p.x) * (mesh.data.pointlist[i].x - p.x) + (mesh.data.pointlist[i].y - p.y) * (mesh.data.pointlist[i].y - p.y);
            if (distance < closest_distance) {
                closest_node = i;
                closest_distance = distance;
            }
        }

        if (closest_node == -1) {
            nlogerror("could not find closest node");
            return;
        }

        closest_node += offsets[layer_id]; // add the offset

        // set the voltage
        system.addDirichletBoundaryCondition(closest_node, V);
    }

    void MultiLayerCurrentDensitySimulation::solve(CurrentDensitySystem& system) {
        V = CV(system.b.val.size());

        MatCSRSymmetric FemMat(system.A);

        auto start = std::chrono::high_resolution_clock::now();
        preconditionedSSORConjugateGradientSolver(FemMat, system.b, V, 1.5, 1e-12, 100000);
        // preconditionedJacobiConjugateGradientSolver(FemMat, system.b, V, 1e-6, 100000);
        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    }

    void MultiLayerCurrentDensitySimulation::Vplot(uint32_t width, uint32_t height) {
        std::vector<cv::Mat> images;
        for (uint64_t i = 0; i < meshes.size(); i++) {
            images.push_back(cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255)));
        }
        
        for (uint64_t i = 0; i < meshes.size(); i++) {
            VplotRend(&images[i], width, height, i);
            cv::flip(images[i], images[i], 0); // flip image horizontally
        }

        // stack the images horizontally and show them
        cv::Mat image;
        cv::hconcat(images, image);
        cv::imshow("V", image);
        cv::waitKey(0);
    }

    void MultiLayerCurrentDensitySimulation::VplotRend(cv::Mat* image, double width, double height, uint64_t layer_id) {
        auto& mesh = meshes[layer_id];
        float min_x = mesh.data.pointlist[0].x;
        float min_y = mesh.data.pointlist[0].y;
        float max_x = mesh.data.pointlist[0].x;
        float max_y = mesh.data.pointlist[0].y;
        for (uint32_t i = 1; i < mesh.data.numberofpoints; i++) {
            if (mesh.data.pointlist[i].x < min_x) {
                min_x = mesh.data.pointlist[i].x;
            }
            if (mesh.data.pointlist[i].y < min_y) {
                min_y = mesh.data.pointlist[i].y;
            }
            if (mesh.data.pointlist[i].x > max_x) {
                max_x = mesh.data.pointlist[i].x;
            }
            if (mesh.data.pointlist[i].y > max_y) {
                max_y = mesh.data.pointlist[i].y;
            }
        }

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * width / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * height / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * width - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * height - 0.5 * (max_y + min_y) * y_scale;

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& mesh : meshes) {
            offsets.push_back(offset);
            offset += mesh.data.numberofpoints;
        }

        // get the range of values for the given layer
        std::vector<float> _V = std::vector<float>(V.val.begin() + offsets[layer_id], V.val.begin() + offsets[layer_id] + meshes[layer_id].data.numberofpoints);

        std::vector<float> V_sorted(V.val.size());
        std::copy(V.val.begin(), V.val.end(), V_sorted.begin());
        std::sort(V_sorted.begin(), V_sorted.end());
        float max_V = V_sorted[0.9 * V_sorted.size()];
        float min_V = V_sorted[0.1 * V_sorted.size()];

        nloginfo("max V: %f", max_V);
        nloginfo("min V: %f", min_V);

        // draw the points
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector p = mesh.data.pointlist[i];

            auto points = std::vector<cv::Point>();

            cv::Scalar c = val2jet(_V[i], min_V, max_V);

            for (int32_t j = 0; j < mesh.data.numberoftriangles; j++) {
                if (mesh.data.trianglelist[j][0] == i || mesh.data.trianglelist[j][1] == i || mesh.data.trianglelist[j][2] == i) {
                    Vector barycenter = {
                        (mesh.data.pointlist[mesh.data.trianglelist[j][0]].x + mesh.data.pointlist[mesh.data.trianglelist[j][1]].x + mesh.data.pointlist[mesh.data.trianglelist[j][2]].x) / 3,
                        (mesh.data.pointlist[mesh.data.trianglelist[j][0]].y + mesh.data.pointlist[mesh.data.trianglelist[j][1]].y + mesh.data.pointlist[mesh.data.trianglelist[j][2]].y) / 3
                    };
                    points.push_back(cv::Point(x_offset + barycenter.x * x_scale, y_offset + barycenter.y * y_scale));
                }
            }

            // find the center of the points
            cv::Point center = {0, 0};
            for (uint8_t j = 0; j < points.size(); j++) {
                center.x += points[j].x;
                center.y += points[j].y;
            }
            // printf("count: %d\n", points.size());
            center.x /= points.size();
            center.y /= points.size();
            // printf("center: %d %d\n", center.x, center.y);

            // sort the points by angle
            std::sort(points.begin(), points.end(), [center](cv::Point a, cv::Point b) {
                return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
            });

            // draw the polygon
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, points), c);
        }
    }
}