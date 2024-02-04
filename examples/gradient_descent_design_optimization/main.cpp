#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

uint64_t time_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

double cost_function(double winding_thickness, double height, double iron_thickness, bool plot = false) {
    auto sim = nikfemm::MagnetostaticSimulation(1, 0.001);

    double winding_radius = (iron_thickness + winding_thickness) / 2;
    double iron_radius = iron_thickness / 2;
    double total_current = 1000;
    double winding_area = winding_thickness * height * 2;
    float current_density = total_current / winding_area;

    sim.mesh.drawing.drawRectangle({iron_radius, 0}, {winding_radius, height});
    sim.mesh.drawing.drawRectangle({-iron_radius, 0}, {-winding_radius, height});
    sim.mesh.drawing.drawRectangle({iron_radius, 0}, {-iron_radius, height});
    sim.mesh.drawing.drawRectangle({-0.5, -1}, {0.5, -2});

    sim.mesh.drawing.drawRegion({0, height / 2}, {0, {0, 0}, nikfemm::magnetostatic_materials::iron_linear});
    sim.mesh.drawing.drawRegion({iron_radius + winding_thickness / 4, height / 2}, {current_density, {0, 0}, nikfemm::magnetostatic_materials::copper});
    sim.mesh.drawing.drawRegion({-iron_radius - winding_thickness / 4, height / 2}, {-current_density, {0, 0}, nikfemm::magnetostatic_materials::copper});
    sim.mesh.drawing.drawRegion({0, -1.5}, {0, {0, 0}, nikfemm::magnetostatic_materials::iron_linear});

    auto system = sim.generateSystem();
    sim.solve(system);

    auto force = sim.computeForceIntegrals({0, -1.5});
#ifdef NIKFEMM_USE_OPENCV
    if (plot) sim.Bplot(1000, 1000, false, true, NAN, NAN, false);
#endif
    return force.y;
}

int main(int argc, char** argv) {
    double epsilon = 0.3;
    double step_size = 1000;
    double winding_thickness = 2;
    double height = 2;
    double iron_thickness = 2;
    // double iron_thickness = 1.9999991421788814 + 0.29999999999999999;
    // double winding_thickness = 1.9999997851909908;
    // double height = 1.9999996387230858;
    for (int i = 0; i < 1000; i++) {
        // compute gradient
        double base = cost_function(winding_thickness, height, iron_thickness, true);
        double dF1 = cost_function(winding_thickness + epsilon, height, iron_thickness) - base;
        double dF2 = cost_function(winding_thickness, height + epsilon, iron_thickness) - base;
        double dF3 = cost_function(winding_thickness, height, iron_thickness + epsilon) - base;

        double dF1d = dF1 / epsilon;
        double dF2d = dF2 / epsilon;
        double dF3d = dF3 / epsilon;

        printf("------------------------------------------------------------------\n");
        printf("Iteration %d\n", i);
        printf("cost: %.5g, winding_thickness: %f, height: %f, iron_thickness: %f, dF1d: %.5g, dF2d: %.5g, dF3d: %.5g\n", base, winding_thickness, height, iron_thickness, dF1d, dF2d, dF3d);

        // update parameters
        winding_thickness += step_size * dF1d;
        // height += step_size * dF2d;
        iron_thickness += step_size * dF3d;
    }
}