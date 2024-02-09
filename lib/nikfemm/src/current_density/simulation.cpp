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

#include "../constants.hpp"
#include "simulation.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/vector.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    CurrentDensitySimulation::CurrentDensitySimulation(double depth) {
        mesh.depth = depth;
    }

    CurrentDensitySimulation::CurrentDensitySimulation() {
        mesh.depth = 1;
        mesh = CurrentDensityMesh();
    }

    CurrentDensitySimulation::~CurrentDensitySimulation() {
        
    }

    void CurrentDensitySimulation::setVoltage(System<double>& system, Vector p, double V) {
        // find the closest node
        int32_t closest_node = -1;
        double closest_distance = INFINITY;
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

        // set the voltage
        system.addDirichletBoundaryCondition(closest_node, V);
    }

    System<double> CurrentDensitySimulation::generateSystem(bool refine, double max_triangle_area, int min_angle) {
        if (refine) {
            mesh.drawing.addRefiningPoints();
        }

        mesh.mesh(max_triangle_area, min_angle);
        mesh.computeEpsilon();
        nloginfo("the mesh has %u nodes and %u elements", mesh.data.numberofpoints, mesh.data.numberoftriangles);

        auto system = mesh.getFemSystem();
        // auto system = mesh.getFemSystemCotangentWeights();

        return system;
    }

    std::vector<double> CurrentDensitySimulation::solve(System<double>& system) {
        auto V = std::vector<double>(mesh.data.numberofpoints);

        MatCSRSymmetric FemMat(system.A);

        auto start = std::chrono::high_resolution_clock::now();
        preconditionedSSORConjugateGradientSolver(FemMat, system.b, V, 1.5, 1e-12, 100000);
        // preconditionedJacobiConjugateGradientSolver(FemMat, system.b, V, 1e-6, 100000);
        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);

        return V;
    }
}