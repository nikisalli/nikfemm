#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::CurrentDensitySimulation simulation;
    simulation.depth = 70e-6;
    simulation.mesh.max_triangle_area = 1e-2;

    // simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(1, 1));

    // simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {nikfemm::current_density_materials::copper});

    simulation.mesh.drawing.drawPolygon(
        {
            nikfemm::Vector(0, 0),
            nikfemm::Vector(4, 0),
            nikfemm::Vector(4, 4),
            nikfemm::Vector(0, 4),
            nikfemm::Vector(0, 2),
            nikfemm::Vector(1, 2),
            nikfemm::Vector(1, 3),
            nikfemm::Vector(3, 3),
            nikfemm::Vector(3, 1),
            nikfemm::Vector(1, 1),
            nikfemm::Vector(1, 1.5),
            nikfemm::Vector(0, 1.5)

        }
    );

    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {nikfemm::current_density_materials::copper});

    auto system = simulation.generateSystem();

    simulation.setVoltage(system, nikfemm::Vector(0.5, 0.5), 1);
    simulation.setVoltage(system, nikfemm::Vector(3.5, 3.5), 0);

    simulation.solve(system);

    simulation.Vplot(1000, 1000);
}