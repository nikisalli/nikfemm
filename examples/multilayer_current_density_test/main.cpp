#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MultiLayerCurrentDensitySimulation simulation(2);
    simulation.depths = {70e-6, 70e-6};
    simulation.meshes[0].max_triangle_area = 1e-3;
    simulation.meshes[1].max_triangle_area = 1e-3;

    // simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(1, 1));

    // simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {nikfemm::current_density_materials::copper});

    simulation.meshes[0].drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(2, 2));
    simulation.meshes[1].drawing.drawRectangle(nikfemm::Vector(1, 1), nikfemm::Vector(3, 3));

    simulation.meshes[0].drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {nikfemm::current_density_materials::copper});
    simulation.meshes[1].drawing.drawRegion(nikfemm::Vector(2.5, 2.5), {nikfemm::current_density_materials::copper});

    simulation.interconnections.push_back({nikfemm::Vector(1.5, 1.5), nikfemm::Vector(1.5, 1.5), 0, 1, 100});

    auto system = simulation.generateSystem();

    simulation.setVoltage(system, nikfemm::Vector(0.1, 0.1), -1, 0);
    simulation.setVoltage(system, nikfemm::Vector(2.9, 2.9), 1, 1);

    simulation.solve(system);

    simulation.Vplot(400, 400);
}