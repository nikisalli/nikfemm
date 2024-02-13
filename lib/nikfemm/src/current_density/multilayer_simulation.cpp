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

#include "../triangle/triangle.h"

#include "../constants.hpp"
#include "multilayer_simulation.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/vector.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths) {
        if (depths.size() != num_layers) {
            throw std::invalid_argument("depths.size() != num_layers");
        }

        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].depth = depths[i];
        }
    }

    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers) {
        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].depth = 1;
        }
    }

    static inline double geomAngle(Vector a, Vector b, Vector c) {
        // safe angle calculation
        double ax = a.x - b.x;
        double ay = a.y - b.y;
        double bx = c.x - b.x;
        double by = c.y - b.y;
        double dot = ax * bx + ay * by;
        double det = ax * by - ay * bx;
        double angle = fabs(atan2(det, dot));
        return angle;
    }

    System<double> MultiLayerCurrentDensitySimulation::generateSystem(bool refine, double max_triangle_area, int min_angle) {
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

        // refine the mesh
        for (auto& interconnection : interconnections) {
            // first we need to compute an adequate epsilon for the ring around the interconnection
            // to compute it we will use the distance from the nearest segment in the drawing
            double epsilon1 = std::numeric_limits<double>::max();
            double epsilon2 = std::numeric_limits<double>::max();
            auto& drawing1 = meshes[interconnection.layer1_id].drawing;
            auto& drawing2 = meshes[interconnection.layer2_id].drawing;

            for (uint64_t i = 0; i < drawing1.segments.size(); i++) {
                double dist = Segment::pointSegmentDistance(interconnection.p1, drawing1.points[drawing1.segments[i].p1], drawing1.points[drawing1.segments[i].p2]);
                if (dist < epsilon1) {
                    epsilon1 = dist;
                }
            }

            for (uint64_t i = 0; i < drawing2.segments.size(); i++) {
                double dist = Segment::pointSegmentDistance(interconnection.p2, drawing2.points[drawing2.segments[i].p1], drawing2.points[drawing2.segments[i].p2]);
                if (dist < epsilon2) {
                    epsilon2 = dist;
                }
            }

            // add points to the mesh all around the interconnection
            const double num_points = 10;
            for (uint32_t i = 0; i < num_points; i++) {
                double angle = 2 * M_PI * i / num_points;
                Vector p1 = interconnection.p1 + Vector(epsilon1 * cos(angle) * 0.5, epsilon1 * sin(angle) * 0.5);
                Vector p2 = interconnection.p2 + Vector(epsilon2 * cos(angle) * 0.5, epsilon2 * sin(angle) * 0.5);
                Vector p1inner = interconnection.p1 + Vector(epsilon1 * cos(angle) * 0.25, epsilon1 * sin(angle) * 0.25);
                Vector p2inner = interconnection.p2 + Vector(epsilon2 * cos(angle) * 0.25, epsilon2 * sin(angle) * 0.25);
                meshes[interconnection.layer1_id].drawing.drawPoint(p1);
                meshes[interconnection.layer2_id].drawing.drawPoint(p2);
                meshes[interconnection.layer1_id].drawing.drawPoint(p1inner);
                meshes[interconnection.layer2_id].drawing.drawPoint(p2inner);
            }
        }

        for (uint64_t i = 0; i < meshes.size(); i++) {
            meshes[i].mesh(max_triangle_area, min_angle);
            meshes[i].computeEpsilon();
            nloginfo("Meshed layer %d with %d nodes and %d elements", i, meshes[i].data.numberofpoints, meshes[i].data.numberoftriangles);
        }

        std::vector<System<double>> systems;

        for (uint64_t i = 0; i < meshes.size(); i++) {
            systems.push_back(meshes[i].getFemSystem());
        }

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& system : systems) {
            offsets.push_back(offset);
            offset += system.b.size();
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

            // nloginfo("Nearest node to p1 is %d at x: %f, y: %f", node1_id, meshes[interconnection.layer1_id].data.pointlist[node1_id].x, meshes[interconnection.layer1_id].data.pointlist[node1_id].y);
            // nloginfo("Nearest node to p2 is %d at x: %f, y: %f", node2_id, meshes[interconnection.layer2_id].data.pointlist[node2_id].x, meshes[interconnection.layer2_id].data.pointlist[node2_id].y);

            interconnection_indices_list.push_back({interconnection.layer1_id, interconnection.layer2_id, node1_id, node2_id, interconnection.R});
        }

        // preallocate the size of the merged system by summing the sizes of the systems
        uint64_t merged_system_nnz = 0;
        uint64_t merged_system_size = 0;
        for (auto& system : systems) {
            merged_system_nnz += system.A.elems.size();
            merged_system_size += system.b.size();
        }

        System<double> merged_system(merged_system_size);

        nloginfo("merged_system_nnz: %d", merged_system_nnz);
        nloginfo("merged_system_size: %d", merged_system_size);
        merged_system.A.elems.reserve(merged_system_nnz);

        // merge the systems together like a block matrix
        offset = 0; // reusing offset variable

        for (uint64_t i = 0; i < meshes.size(); i++) {
            System<double>& system = systems[i];
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
                merged_system.b[j + offset] = system.b[j];
            }

            offset += system_size;
        }

        nloginfo("b.size(): %d", merged_system.b.size());

        // add the interconnections to the merged system
        // method 1 (couldn't bother to get the math right)
        for (auto& interconnection : interconnection_indices_list) {
            // add the interconnection to the merged system
            uint64_t m1 = interconnection.node1_id + offsets[interconnection.mesh1_id];
            uint64_t m2 = interconnection.node2_id + offsets[interconnection.mesh2_id];
            double R = interconnection.R;

            // 1.
            merged_system.A(m1, m1) += 2/R;
            merged_system.A(m2, m2) += 2/R;

            // 2.
            if (m1 >= m2) {
                merged_system.A(m2, m1) = -2/R;
            } else {
                merged_system.A(m1, m2) = -2/R;
            }
        }

        return merged_system;
    }

    void MultiLayerCurrentDensitySimulation::setVoltage(System<double>& system, Vector p, double V, uint64_t layer_id) {
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

    std::vector<double> MultiLayerCurrentDensitySimulation::solve(System<double>& system) {
        auto V = std::vector<double>(system.b.size());

        MatCSRSymmetric FemMat(system.A);

        auto start = std::chrono::high_resolution_clock::now();
        preconditionedSSORConjugateGradientSolver(FemMat, system.b, V, 1.5, 1e-7, 100000);
        // preconditionedJacobiConjugateGradientSolver(FemMat, system.b, V, 1e-6, 100000);
        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);

        return V;
    }
}