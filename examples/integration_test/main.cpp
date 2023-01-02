#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation;
   
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(1, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, -2), nikfemm::Vector(1, -1));

    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {0, {0, 100}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, -1.5), {0, {0, -100}, nikfemm::materials::air});

    auto system = simulation.generateSystem();
    simulation.solve(system);
    // simulation.Aplot(400, 400);
    // simulation.Bplot(400, 400);
    // simulation.AplotToFile(100000, 100000, "Aplot.png");
    simulation.BplotToFile(10000, 10000, "Bplot.png", false, false);

    // auto force = simulation.computeForceIntegrals({0.5, 0.5});
    // printf("Force: %.17g, %.17g\n", force.x, force.y);
}