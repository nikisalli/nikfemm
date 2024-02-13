#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation;

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
    
   
    /*
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(0, 0), nikfemm::Vector(1, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(1, 0), nikfemm::Vector(2, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(2, 0), nikfemm::Vector(3, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(3, 0), nikfemm::Vector(4, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Vector(4, 0), nikfemm::Vector(5, 1));

    simulation.mesh.drawing.drawRegion(nikfemm::Vector(-1, 0), {0, {0, 0}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(0.5, 0.5), {0, {-1, 0}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(1.5, 0.5), {0, {0, -1}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(2.5, 0.5), {0, {1, 0}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(3.5, 0.5), {0, {0, 1}, nikfemm::magnetostatic_materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Vector(4.5, 0.5), {0, {-1, 0}, nikfemm::magnetostatic_materials::air});
    */

    auto system = simulation.generateSystem();
    auto A = simulation.solve(system);
    auto B = simulation.mesh.computeCurl(A);
    // simulation.Aplot(400, 400);
    // simulation.Bplot(400, 400);
    // simulation.AplotToFile(100000, 100000, "Aplot.png");
    // simulation.BplotToFile(10000, 10000, "Bplot.png", false, false);
#ifdef NIKFEMM_USE_OPENCV
    // simulation.mesh.ElemScalarPlot(1000, 1000, B, false, false, true);
#endif

    // simulation.computeForceIntegrals({0.5, 0.5});
}