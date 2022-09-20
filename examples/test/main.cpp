#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::Simulation simulation;

    nikfemm::Drawing drawing;

    drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 1));
    // drawing.drawRectangle(nikfemm::Point(0, 1), nikfemm::Point(1, 2));
    drawing.drawRegion(nikfemm::Point(0.5, 0.5), 6000);

    simulation.generateMesh(drawing);
}