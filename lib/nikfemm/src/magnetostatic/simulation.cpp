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
#include "../geometry/point.hpp"
#include "../algebra/coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    MagnetostaticSimulation::MagnetostaticSimulation() {

    }

    MagnetostaticSimulation::~MagnetostaticSimulation() {

    }

    void MagnetostaticSimulation::generateMesh() {
        // get time in milliseconds
        mesh.drawing.addRefiningPoints();

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        auto start1 = std::chrono::high_resolution_clock::now();
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(mesh.drawing.points);
        if (mesh.drawing.points.size() == 0) {
            smallest_circle.radius = 1;
        }
        // translate everything to the origin
        for (uint32_t i = 0; i < mesh.drawing.points.size(); i++) {
            mesh.drawing.points[i] = Point(mesh.drawing.points[i].x - smallest_circle.center.x, mesh.drawing.points[i].y - smallest_circle.center.y);
        }
        std::vector<DrawingRegion> translated_regions;
        for (DrawingRegion region : mesh.drawing.regions) {
            translated_regions.push_back(DrawingRegion(Point(region.first.x - smallest_circle.center.x, region.first.y - smallest_circle.center.y), region.second));
        }
        mesh.drawing.regions = translated_regions;
        // set simulation offset and boundary radius
        mesh.center = smallest_circle.center;
        mesh.radius = smallest_circle.radius * 2;
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(Point(0, 0), 2 * smallest_circle.radius);
        mesh.drawing.drawCircle(boundary_circle, BOUNDARY_VERTICES);
        // add region near the edge of the circle
        mesh.drawing.drawRegion(Point(boundary_circle.radius * 0.9, 0), {0, {0, 0}, materials::air});
        // add the boundary 
        // mesh.drawing.plot();
        auto start2 = std::chrono::high_resolution_clock::now();
        mesh.mesh();
        auto start3 = std::chrono::high_resolution_clock::now();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        mesh.addKelvinBoundaryConditions();
        auto start4 = std::chrono::high_resolution_clock::now();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        printf("the mesh has %lu nodes and %lu elements\n", mesh.data.numberofpoints, mesh.data.numberoftriangles);
        MatCOO<MagnetostaticNonLinearExpression> coo(mesh.data.numberofpoints);
        CV b(mesh.data.numberofpoints);
        CV x(mesh.data.numberofpoints);
        std::vector<Vector> B(mesh.data.numberoftriangles, {0, 0});
        std::vector<float> mu(mesh.data.numberoftriangles, 0);
        std::vector<const MagnetostaticProp*> props(mesh.data.numberoftriangles);

        // fill props
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            props[i] = mesh.drawing.getRegionPtrFromId(mesh.data.triangleattributelist[i]);
        }

        auto start5 = std::chrono::high_resolution_clock::now();

        mesh.getFemSystem(coo, b);
        auto start6 = std::chrono::high_resolution_clock::now();
        mesh.addDirichletBoundaryConditions(coo, b);
        auto start7 = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        printf("coo matrix m: %lu, n: %lu, elems: %lu\n", coo.m, coo.n, coo.elems.size());
        #endif

        MagnetostaticMatCSRSymmetric A(coo);
        // b.write_to_file("b");
        auto start8 = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        // csr.print();
        // coo.plot();
        // x.print();
        // b.print();
        #endif

        // initialize mu
        for (uint32_t i = 0; i < B.size(); i++) {
            mu[i] = props[i]->getMu(0);
        }
        A.updateMat(mu);

        CV r(b.m);  // residual

        double old_residual = 0;
        double residual = 1e10;
        for (uint32_t i = 0; i < 500; i++) {
            // conjugateGradientSolver(A, b, x, 1e-7, 10000);
            // preconditionedJacobiConjugateGradientSolver(A, b, x, 1e-7, 1000);
            preconditionedSSORConjugateGradientSolver(A, b, x, 1.5, 1e-7, 50);
            // preconditionedIncompleteCholeskyConjugateGradientSolver(A, b, x, 1e-7, 1000);
            // mesh.Aplot(x);
            // mesh.Bplot(B);

            mesh.computeCurl(B, x);
            // mesh.Bplot(B);

            // check if the solution is correct
            A.updateMu(props, mu, B, residual, i);
            A.updateMat(mu);
            CV::mult(r, A, x);
            CV::sub(r, b, r);
            old_residual = residual;
            residual = CV::norm(r);
            if (residual < 1e-7) {
                printf("Converged in %d iterations\n", i);
                break;
            }
            // printf("%.17g, %.17g; ", K, residual);
            printf("nonlinear iteration %d, residual: %.17g\n", i, residual);
            // print mu
            // printf("%f,", residual);
            fflush(stdout);
        }

        auto start9 = std::chrono::high_resolution_clock::now();

        mesh.computeCurl(B, x);

        auto end = std::chrono::high_resolution_clock::now();
        printf("%f translate and fix mesh\n", std::chrono::duration_cast<std::chrono::duration<float>>(start2 - start1).count()*1000);
        printf("%f mesh\n", std::chrono::duration_cast<std::chrono::duration<float>>(start3 - start2).count()*1000);
        printf("%f kelvin boundary conditions\n", std::chrono::duration_cast<std::chrono::duration<float>>(start4 - start3).count()*1000);
        printf("%f allocate vectors b x and coo\n", std::chrono::duration_cast<std::chrono::duration<float>>(start5 - start4).count()*1000);
        printf("%f get fem system\n", std::chrono::duration_cast<std::chrono::duration<float>>(start6 - start5).count()*1000);
        printf("%f add dirichlet boundary conditions\n", std::chrono::duration_cast<std::chrono::duration<float>>(start7 - start6).count()*1000);
        printf("%f convert to csr\n", std::chrono::duration_cast<std::chrono::duration<float>>(start8 - start7).count()*1000);
        printf("%f solve\n", std::chrono::duration_cast<std::chrono::duration<float>>(start9 - start8).count()*1000);
        printf("%f total\n", std::chrono::duration_cast<std::chrono::duration<float>>(start9 - start1).count()*1000);

        mesh.Aplot(x);
        mesh.Bplot(B);
        // mesh.muplot(mu);
        // return;
        /*
        mesh.Bplot();

        // report(&out, 1, 1, 0, 0, 0, 0);

        printf("Number of points: %d\nNumber of boundary vertices: %d\n", mesh.vertices.size(), mesh.boundary_vertices.size());
        */

    }
}