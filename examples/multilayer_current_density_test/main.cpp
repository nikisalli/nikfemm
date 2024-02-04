#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MultiLayerCurrentDensitySimulation simulation(2);
    simulation.meshes[0].depth = 70e-6;
    simulation.meshes[1].depth = 70e-6;
    simulation.meshes[0].max_triangle_area = 1;
    simulation.meshes[1].max_triangle_area = 1;

    // two long thin rectangles coupled together at the ends
    double length = 1;
    double width = 0.001;
    simulation.meshes[0].drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(length, width));
    simulation.meshes[1].drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(length, width));

    simulation.meshes[0].drawing.drawRegion(nikfemm::Vector(0.5 * length, 0.5 * width), {nikfemm::current_density_materials::copper});
    simulation.meshes[1].drawing.drawRegion(nikfemm::Vector(0.5 * length, 0.5 * width), {nikfemm::current_density_materials::copper});

    simulation.interconnections.push_back({nikfemm::Vector(length - 0.5 * width, 0.5 * width), nikfemm::Vector(0.5 * width, 0.5 * width), 0, 1, 1});
    
    auto system = simulation.generateSystem();

    simulation.setVoltage(system, nikfemm::Vector(0, 0.5 * width), -1, 0);
    simulation.setVoltage(system, nikfemm::Vector(length, 0.5 * width), 1, 1);

    simulation.solve(system);

#ifdef NIKFEMM_USE_OPENCV
    simulation.Vplot(1000, 1000);
#endif
}