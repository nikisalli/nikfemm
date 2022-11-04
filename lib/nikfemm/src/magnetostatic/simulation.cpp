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

namespace nikfemm {
    MagnetostaticSimulation::MagnetostaticSimulation() {

    }

    MagnetostaticSimulation::~MagnetostaticSimulation() {

    }

    void MagnetostaticSimulation::generateMesh() {
        // get time in milliseconds

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(mesh.drawing.points);
        if (mesh.drawing.points.size() == 0) {
            smallest_circle.radius = 1;
        }
        // translate everything to the origin
        for (uint64_t i = 0; i < mesh.drawing.points.size(); i++) {
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
        mesh.drawing.drawRegion(Point(boundary_circle.radius * 0.9, 0), vacuum_prop);
        // add the boundary 
        // mesh.drawing.plot();
        auto start = std::chrono::high_resolution_clock::now();
        mesh.mesh();
#ifdef DEBUG_PRINT
        mesh.plot();
#endif
        mesh.addKelvinBoundaryConditions();
#ifdef DEBUG_PRINT
        mesh.plot();
#endif
        MatCOO coo;
        CV b(mesh.data.numberofpoints);
        CV x(mesh.data.numberofpoints);

        // initialize x to 0
        for (uint64_t i = 0; i < mesh.data.numberofpoints; i++) {
            x[i] = 0;
            b[i] = 0;
        }

        mesh.getFemSystem(coo, b);
        mesh.addDirichletBoundaryConditions(coo, b);

#ifdef DEBUG_PRINT
        printf("coo matrix m: %lu, n: %lu, elems: %lu\n", coo.m, coo.n, coo.elems.size());
#endif

        MatCSR csr(coo);
        csr.write_to_file("A");
        b.write_to_file("b");

#ifdef DEBUG_PRINT
        // csr.print();
        coo.plot();
        // x.print();
        // b.print();
#endif

        csr.conjugateGradientSolve(b, x, 1e-7, 10000);

        mesh.Aplot(x);
        return;
        /*
        mesh.Bplot();

        // report(&out, 1, 1, 0, 0, 0, 0);

        printf("Number of points: %d\nNumber of boundary vertices: %d\n", mesh.vertices.size(), mesh.boundary_vertices.size());
        */
        // get end time
        auto end = std::chrono::high_resolution_clock::now();
        printf("%f Total time\n", std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count()*1000);
    }
}