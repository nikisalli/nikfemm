#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <iterator>

#include <simulation.hpp>
#include <drawing_region.hpp>
#include <drawing.hpp>

#include "geometry/segment.hpp"
#include "geometry/point.hpp"

namespace nikfemm {
    Simulation::Simulation() {

    }

    Simulation::~Simulation() {

    }

    void Simulation::generateMesh(Drawing drawing) {
        // get time in milliseconds
        auto start = std::chrono::high_resolution_clock::now();

        std::unordered_set<Point> valid_points;
        for (auto p : drawing.segments) {
            valid_points.insert(p.p1);
            valid_points.insert(p.p2);
        }

        auto start2 = std::chrono::high_resolution_clock::now();

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(valid_points);
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(smallest_circle.center, 2 * smallest_circle.radius);
        drawing.drawCircle(boundary_circle, 360);
        // add region near the edge of the circle
        drawing.drawRegion(Point(boundary_circle.center.x + boundary_circle.radius - EPSILON, boundary_circle.center.y), BOUNDARY_REGION);
        // add the boundary 

        auto start3 = std::chrono::high_resolution_clock::now();

        mesh = drawing.mesh();

        auto start6 = std::chrono::high_resolution_clock::now();

        // report(&out, 1, 1, 0, 0, 0, 0);

        printf("Number of points: %d\nNumber of triangles: %d\nNumber of boundary vertices: %d\nNumber of boundary elements: %d\n", mesh->vertices.size(), mesh->elements.size(), mesh->boundary_vertices.size(), mesh->boundary_elements.size());
        // get end time
        auto end = std::chrono::high_resolution_clock::now();
        printf("%f Times:\n");
        printf("%f Drawing segments to Verts and Segs\n", std::chrono::duration_cast<std::chrono::duration<double>>(start2 - start).count()*1000);
        printf("%f Draw boundary\n", std::chrono::duration_cast<std::chrono::duration<double>>(start3 - start2).count()*1000);
        printf("%f Triangulate\n", std::chrono::duration_cast<std::chrono::duration<double>>(start6 - start3).count()*1000);
        printf("%f Find boundary vertices\n", std::chrono::duration_cast<std::chrono::duration<double>>(end - start6).count()*1000);
        printf("%f Total time\n", std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count()*1000);
        
        // mesh->plot();
    }
}