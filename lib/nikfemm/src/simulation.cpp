#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <iterator>

#include <constants.hpp>
#include <simulation.hpp>

#include "drawing/drawing_region.hpp"
#include "drawing/drawing.hpp"
#include "geometry/segment.hpp"
#include "geometry/point.hpp"

#include "matrix/coo.hpp"
#include "matrix/csr.hpp"
#include "matrix/simple_vector.hpp"

namespace nikfemm {
    Simulation::Simulation() {

    }

    Simulation::~Simulation() {

    }

    void Simulation::generateMesh(Drawing drawing) {
        // get time in milliseconds

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(drawing.points);
        // translate everything to the origin
        for (uint64_t i = 0; i < drawing.points.size(); i++) {
            drawing.points[i] = Point(drawing.points[i].x - smallest_circle.center.x, drawing.points[i].y - smallest_circle.center.y);
        }
        std::unordered_set<DrawingRegion> translated_regions;
        for (DrawingRegion region : drawing.regions) {
            translated_regions.insert(DrawingRegion(Point(region.p.x - smallest_circle.center.x, region.p.y - smallest_circle.center.y), region.region_attribute));
        }
        drawing.regions = translated_regions;
        // set simulation offset and boundary radius
        mesh.center = smallest_circle.center;
        mesh.radius = smallest_circle.radius * 2;
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(Point(0, 0), 2 * smallest_circle.radius);
        drawing.drawCircle(boundary_circle, BOUNDARY_VERTICES);
        // add region near the edge of the circle
        drawing.drawRegion(Point(boundary_circle.radius * 0.9, 0), BOUNDARY_REGION);
        // add the boundary 
        drawing.plot();
        mesh.center = boundary_circle.center;
        auto start = std::chrono::high_resolution_clock::now();
        mesh.mesh(drawing);
        mesh.plot();
        mesh.addKelvinBoundaryConditions();
        mesh.plot();
        MatCOO coo;
        mesh.getFemMatrix(coo);

        // report(&out, 1, 1, 0, 0, 0, 0);

        printf("Number of points: %d\nNumber of boundary vertices: %d\n", mesh.vertices.size(), mesh.boundary_vertices.size());
        // get end time
        auto end = std::chrono::high_resolution_clock::now();
        printf("%f Total time\n", std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count()*1000);
    }
}