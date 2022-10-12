#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation;

    // drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 1));
    // drawing.drawRectangle(nikfemm::Point(0, 1), nikfemm::Point(1, 2));

    // drawing.drawCircle(nikfemm::Point(1, 0), 0.5, 100);
    // drawing.drawRegion(nikfemm::Point(1, 0), 1);
    // drawing.drawCircle(nikfemm::Point(-1, 0), 0.5, 100);
    // drawing.drawRegion(nikfemm::Point(-1, 0), -1);

    simulation.drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 2));
    simulation.drawing.drawRegion(nikfemm::Point(0.5, 0.5), {MU_0, -1, {0, 0}});
    simulation.drawing.drawRectangle(nikfemm::Point(0, 3), nikfemm::Point(1, 5));
    simulation.drawing.drawRegion(nikfemm::Point(0.5, 4), {MU_0, -1, {0, 0}});
    simulation.drawing.drawRectangle(nikfemm::Point(4, 0), nikfemm::Point(5, 5));
    simulation.drawing.drawRegion(nikfemm::Point(4.5, 2.5), {MU_0, 1, {0, 0}});

    simulation.generateMesh();
}