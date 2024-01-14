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
            }

            for (uint64_t i = 0; i < meshes[interconnection.layer2_id].data.numberofpoints; i++) {
                double dist = 0;

                dist = (meshes[interconnection.layer2_id].data.pointlist[i] - interconnection.p2).norm();
                if (dist < node2_dist) {
                    node2_id = i;
                    node2_dist = dist;
                }
            }

            nloginfo("Nearest node to p1 is %d at x: %f, y: %f", node1_id, meshes[interconnection.layer1_id].data.pointlist[node1_id].x, meshes[interconnection.layer1_id].data.pointlist[node1_id].y);
            nloginfo("Nearest node to p2 is %d at x: %f, y: %f", node2_id, meshes[interconnection.layer2_id].data.pointlist[node2_id].x, meshes[interconnection.layer2_id].data.pointlist[node2_id].y);

            interconnection_indices_list.push_back({interconnection.layer1_id, interconnection.layer2_id, node1_id, node2_id, interconnection.R});
        }

        // preallocate the size of the merged system by summing the sizes of the systems
        uint64_t merged_system_nnz = 0;
        uint64_t merged_system_size = 0;
        nloginfo("merged_system_nnz is: %d", merged_system_nnz);
        nloginfo("merged_system_size is: %d", merged_system_size);
        for (auto& system : systems) {
            merged_system_nnz += system.A.elems.size();
            merged_system_size += system.b.val.size();
        }

        CurrentDensitySystem merged_system = {
            BuildMatCOO<double>(merged_system_size),
            CV(merged_system_size)
        };

        nloginfo("merged_system_nnz: %d", merged_system_nnz);
        nloginfo("merged_system_size: %d", merged_system_size);
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

                merged_system.A(m + offset, n + offset) = elem.second;
            }

            // b
            for (uint64_t j = 0; j < system_size; j++) {
                merged_system.b.val.push_back(system.b.val[j]);
            }

            offset += system_size;
        }

        nloginfo("b.size(): %d", merged_system.b.val.size());

        // add the interconnections to the merged system        
        for (auto& interconnection : interconnection_indices_list) {
            // add the interconnection to the merged system
            uint64_t m1 = interconnection.node1_id + offsets[interconnection.mesh1_id];
            uint64_t m2 = interconnection.node2_id + offsets[interconnection.mesh2_id];
            double R = interconnection.R;

            // get all the rows that have a nonzero element in column m1
            std::unordered_set<uint64_t> rows_with_m1;
            std::unordered_set<uint64_t> rows_with_m2;

            for (auto& elem : merged_system.A.elems) {
                uint64_t m = elem.first >> 32;
                uint64_t n = elem.first & 0xFFFFFFFF;

                if (n == m1) {
                    rows_with_m1.insert(m);
                }
                if (n == m2) {
                    rows_with_m2.insert(m);
                }
            }

            // now for every row that has a nonzero element in column m1 we should:
            // 1. for every nonzero element in that row, multiply it by R and add + 1 if the column index is m1
            // 2. add a new element to the merged system with the value -1 in column m2 and the same row index

            // then for every row that has a nonzero element in column m2 we should:
            // 1. for every nonzero element in that row, multiply it by R and add + 1 if the column index is m2
            // 2. add a new element to the merged system with the value -1 in column m1 and the same row index

            // 1.
            for (auto& elem : merged_system.A.elems) {
                uint64_t m = elem.first >> 32;
                uint64_t n = elem.first & 0xFFFFFFFF;

                if (rows_with_m1.find(m) != rows_with_m1.end()) {
                    elem.second *= (-1 / R) + 1;
                    if (n == m1) {
                        elem.second += 1;
                    }
                } else if (rows_with_m2.find(m) != rows_with_m2.end()) {
                    elem.second *= (-1 / R) + 1;
                    if (n == m2) {
                        elem.second += 1;
                    }
                }
            }

            // 2.
            for (auto& elem : merged_system.A.elems) {
                uint64_t m = elem.first >> 32;
                uint64_t n = elem.first & 0xFFFFFFFF;

                if (rows_with_m1.find(m) != rows_with_m1.end()) {
                    if (n == m1) {
                        // the system is symmetric so we only need to add one of the two elements
                        // since it is stored in upper triangular form
                        if(m2 >= m) merged_system.A(m, m2) = +1;
                        else merged_system.A(m2, m) = +1;
                    }
                } else if (rows_with_m2.find(m) != rows_with_m2.end()) {
                    if (n == m2) {
                        // the system is symmetric so we only need to add one of the two elements
                        // since it is stored in upper triangular form
                        if(m1 >= m) merged_system.A(m, m1) = +1;
                        else merged_system.A(m1, m) = +1;
                    }
                }
            }
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
        for (auto& _mesh : meshes) {
            offsets.push_back(offset);
            offset += _mesh.data.numberofpoints;
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
        preconditionedSSORConjugateGradientSolver(FemMat, system.b, V, 1.5, 1e-10, 100000);
        // preconditionedJacobiConjugateGradientSolver(FemMat, system.b, V, 1e-6, 100000);
        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    }

    void MultiLayerCurrentDensitySimulation::Vplot(uint32_t width, uint32_t height) {
        std::vector<cv::Mat> images;
        for (uint64_t i = 0; i < meshes.size(); i++) {
            images.push_back(cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255)));
        }

        // get 1st and 99th percentile of V
        std::vector<double> _V = std::vector<double>(V.val.begin(), V.val.end());
        std::sort(_V.begin(), _V.end());
        double max_V = _V[_V.size() * 0.9];
        double min_V = _V[_V.size() * 0.1];

        max_V = 1;
        min_V = -1;

        nloginfo("max_V: %f", max_V);
        nloginfo("min_V: %f", min_V);

        // print middle left and right
        nloginfo("V at middle left: %f", _V[_V.size() * 0.3]);
        nloginfo("V at middle right: %f", _V[_V.size() * 0.7]);

        for (uint64_t i = 0; i < meshes.size(); i++) {
            VplotRend(&images[i], width, height, i, max_V, min_V);
            cv::flip(images[i], images[i], 0); // flip image horizontally
        }

        // stack the images horizontally and show them
        cv::Mat image;
        cv::hconcat(images, image);
        cv::imshow("V", image);
        cv::waitKey(0);
    }

    void MultiLayerCurrentDensitySimulation::VplotRend(cv::Mat* image, double width, double height, uint64_t layer_id, double max_V, double min_V) {
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
            center.x /= points.size();
            center.y /= points.size();

            // sort the points by angle
            std::sort(points.begin(), points.end(), [center](cv::Point a, cv::Point b) {
                return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
            });

            // draw the polygon
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, points), c);
        }
    }
}