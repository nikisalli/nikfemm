#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::CurrentDensitySimulation simulation;
    simulation.mesh.depth = 70e-6;

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

    auto system = simulation.generateSystem(false, 1e-3);

    simulation.setVoltage(system, nikfemm::Vector(0.5, 0.5), 1);
    simulation.setVoltage(system, nikfemm::Vector(3.5, 3.5), 0);

    auto V = simulation.solve(system);

#ifdef NIKFEMM_USE_OPENCV
    simulation.mesh.NodeScalarPlot(1000, 1000, V, false, false, false);
#endif
}