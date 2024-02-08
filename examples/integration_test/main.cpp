#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation;
    simulation.depth = 1.0;
   
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, -0.1), nikfemm::Vector(1, 0.9));
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, 2.1), nikfemm::Vector(1, 3.1));
    simulation.mesh.drawing.drawPolygon(
        {
            nikfemm::Vector(-2, 1),
            nikfemm::Vector(3, 1),
            nikfemm::Vector(3, 4),
            nikfemm::Vector(2, 4),
            nikfemm::Vector(2, 2),
            nikfemm::Vector(-1, 2),
            nikfemm::Vector(-1, 4),
            nikfemm::Vector(-2, 4)
        }
    );

    simulation.mesh.drawing.drawPolygon(
        {
            nikfemm::Vector(3, 4.5),
            nikfemm::Vector(3, 5.5),
            nikfemm::Vector(-2, 5.5),
            nikfemm::Vector(-2, 4.5)
        }
    );

    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {-1, {0, 0}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 2.5), {1, {0, 0}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0, 1.5), {0, {0, 0}, nikfemm::magnetostatic_materials::iron_linear});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0, 5), {0, {0, 0}, nikfemm::magnetostatic_materials::iron_linear});

    auto system = simulation.generateSystem(true, 0.1);
    auto A = simulation.solve(system);
    auto B = simulation.mesh.computeCurl(A);
    // simulation.Aplot(400, 400);
#ifdef NIKFEMM_USE_OPENCV
    simulation.mesh.ElemScalarPlot(1000, 1000, B, false, false, true);
    // simulation.AplotToFile(100000, 100000, "Aplot.png");
    simulation.mesh.ElemScalarPlotToFile(1000, 1000, B, "Bplot.png", false, false);
#endif

    auto stress = simulation.computeStressIntegral(B, {0, 1.5}, {0, 1.5});
    printf("Force: %.17g, %.17g, Torque: %.17g\n", stress.Force.x, stress.Force.y, stress.Torque);
}