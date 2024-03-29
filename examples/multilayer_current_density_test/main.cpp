#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MultiLayerCurrentDensitySimulation simulation(2);
    simulation.meshes[0].depth = 70e-6;
    simulation.meshes[1].depth = 70e-6;

    // two long thin rectangles coupled together at the ends
    double length = 1;
    double width = 0.1;
    simulation.meshes[0].drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(length, width));
    simulation.meshes[1].drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(length, width));

    simulation.meshes[0].drawing.drawRegion(nikfemm::Vector(0.5 * length, 0.5 * width), {nikfemm::current_density_materials::copper});
    simulation.meshes[1].drawing.drawRegion(nikfemm::Vector(0.5 * length, 0.5 * width), {nikfemm::current_density_materials::copper});

    simulation.interconnections.push_back({nikfemm::Vector(length - 0.5 * width, 0.5 * width), nikfemm::Vector(0.5 * width, 0.5 * width), 0, 1, 1});
    
    nikfemm::System<double> system = simulation.generateSystem(false, 0.000001, 20);

    simulation.setVoltage(system, nikfemm::Vector(0, 0.5 * width), -1, 0);
    simulation.setVoltage(system, nikfemm::Vector(length, 0.5 * width), 1, 1);

    auto V = simulation.solve(system);
    auto P = simulation.computePowerDensity(V);

    simulation.meshes[0].ElemScalarPlot(1000, 1000, P, false, false, false);

#ifdef NIKFEMM_USE_OPENCV
    // simulation.meshes[0].NodeScalarPlot(1000, 1000, V, false, false, false);
#endif
}