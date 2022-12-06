#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation;

    /*
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(0, -0.1), nikfemm::Point(1, 0.9));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(0, 2.1), nikfemm::Point(1, 3.1));
    simulation.mesh.drawing.drawPolygon(
        {
            nikfemm::Point(-2, 1),
            nikfemm::Point(3, 1),
            nikfemm::Point(3, 4),
            nikfemm::Point(2, 4),
            nikfemm::Point(2, 2),
            nikfemm::Point(-1, 2),
            nikfemm::Point(-1, 4),
            nikfemm::Point(-2, 4)
        }
    );

    simulation.mesh.drawing.drawPolygon(
        {
            nikfemm::Point(3, 4.5),
            nikfemm::Point(3, 5.5),
            nikfemm::Point(-2, 5.5),
            nikfemm::Point(-2, 4.5)
        }
    );

    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.5, 0.5), {-1, {0, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.5, 2.5), {1, {0, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0, 1.5), {0, {0, 0}, 0, nikfemm::materials::iron});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0, 5), {0, {0, 0}, 0, nikfemm::materials::iron});
    */
    
   
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(0, 0), nikfemm::Point(1, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(1, 0), nikfemm::Point(2, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(2, 0), nikfemm::Point(3, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(3, 0), nikfemm::Point(4, 1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(4, 0), nikfemm::Point(5, 1));

    simulation.mesh.drawing.drawRegion(nikfemm::Point(-1, 0), {0, {0, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.5, 0.5), {0, {-1, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(1.5, 0.5), {0, {0, -1}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(2.5, 0.5), {0, {1, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(3.5, 0.5), {0, {0, 1}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(4.5, 0.5), {0, {-1, 0}, nikfemm::materials::air});

    /*
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(1, 1), nikfemm::Point(0.9, -1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(-1, 1), nikfemm::Point(-0.9, -1));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(1, 5), nikfemm::Point(0.9, 3));
    simulation.mesh.drawing.drawRectangle(nikfemm::Point(-1, 5), nikfemm::Point(-0.9, 3));
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.95, 0), {1, {0, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(-0.95, 0), {-1, {0, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(0.95, 4), {1, {0, 0}, nikfemm::materials::air});
    simulation.mesh.drawing.drawRegion(nikfemm::Point(-0.95, 4), {-1, {0, 0}, nikfemm::materials::air});
    */

    simulation.solve();
    // simulation.Aplot(400, 400);
    // simulation.Bplot(400, 400);
    // simulation.AplotToFile(100000, 100000, "Aplot.png");
    simulation.BplotToFile(10000, 10000, "Bplot.png");
}