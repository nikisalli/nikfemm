#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::CurrentDensitySimulation simulation;
    simulation.depth = 70e-6;
    simulation.mesh.max_triangle_area = 1e-4;

    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(1, 1));

    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {nikfemm::current_density_materials::copper});

    auto system = simulation.generateSystem();

    system.addDirichletBoundaryCondition(0, 0.0);
    system.addDirichletBoundaryCondition(1, 1.0);

    simulation.solve(system);

    simulation.Vplot(1000, 1000);
}