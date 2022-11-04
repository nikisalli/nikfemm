#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation;

    // simulation.drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 1));
    // simulation.drawing.drawRectangle(nikfemm::Point(0, 1), nikfemm::Point(1, 2));

    // simulation.mesh.drawing.drawCircle(nikfemm::Point(0, 0), 1, 100);
    // simulation.mesh.drawing.drawRegion(nikfemm::Point(0, 0), {MU_0, 1, {0, 0}});

    simulation.mesh.drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(0, 2), nikfemm::Point(1, 3));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(-3, 1), nikfemm::Point(3, 2));
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.5, 0.5), {MU_0, 1, {0, 0}});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.5, 2.5), {MU_0, -1, {0, 0}});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0, 1.5), {1, 1, {0, 0}});

    // simulation.drawing.drawCircle(nikfemm::Point(-1, 0), 0.5, 100);
    // simulation.drawing.drawRegion(nikfemm::Point(-1, 0), -1);

    // simulation.drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 1));
    // simulation.drawing.drawRegion(nikfemm::Point(0.5, 0.5), {1, 1, {0, 0}});

    // simulation.drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 2));
    // simulation.drawing.drawRegion(nikfemm::Point(0.5, 0.5), {MU_0, -1, {0, 0}});
    // simulation.drawing.drawRectangle(nikfemm::Point(0, 3), nikfemm::Point(1, 5));
    // simulation.drawing.drawRegion(nikfemm::Point(0.5, 4), {MU_0, -1, {0, 0}});
    // simulation.drawing.drawRectangle(nikfemm::Point(4, 0), nikfemm::Point(5, 5));
    // simulation.drawing.drawRegion(nikfemm::Point(4.5, 2.5), {MU_0, 1, {0, 0}});

    simulation.generateMesh();
}