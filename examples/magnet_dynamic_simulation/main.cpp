#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    nikfemm::MagnetostaticSimulation simulation(1, 0.1);

    nikfemm::Vector center1 = {0, 0};
    nikfemm::Vector center2 = {2, 2};
    double angle1 = 0;
    double angle2 = 0;
    nikfemm::Vector magnetization1 = {0, 1000};
    nikfemm::Vector magnetization2 = {0, 1000};

    // dummy simulation to get moments and masses
    auto square = nikfemm::Polygon(
        {
            nikfemm::Vector(-0.5, -0.5),
            nikfemm::Vector(0.5, -0.5),
            nikfemm::Vector(0.5, 0.5),
            nikfemm::Vector(-0.5, 0.5)
        }
    );

    simulation.mesh.drawing.drawPolygon(square.copy().translate(center1).rotate(angle1, center1));
    simulation.mesh.drawing.drawPolygon(square.copy().translate(center2).rotate(angle2, center2));
    simulation.mesh.drawing.drawRegion(center1, {0, magnetization1.rotate(angle1), nikfemm::materials::neodymium});
    simulation.mesh.drawing.drawRegion(center2, {0, magnetization2.rotate(angle2), nikfemm::materials::neodymium});

    auto system = simulation.generateSystem();
    simulation.solve(system);

    simulation.Bplot(1000, 1000, false, true);

    double moment_of_inertia1 = simulation.computeInertiaMoment(center1, center1, 7);
    double moment_of_inertia2 = simulation.computeInertiaMoment(center2, center2, 7);
    double mass1 = simulation.computeMass(center1, 7);
    double mass2 = simulation.computeMass(center2, 7);
    printf("Moment of inertia 1: %.17g, Moment of inertia 2: %.17g, Mass 1: %.17g, Mass 2: %.17g\n", moment_of_inertia1, moment_of_inertia2, mass1, mass2);

    double dt = 0.1;

    nikfemm::Vector velocity1 = {0, 0};
    nikfemm::Vector velocity2 = {0, 0};
    double angular_velocity1 = 0;
    double angular_velocity2 = 0;
    for (int i = 0; i < 10000; i++) {
        nikfemm::MagnetostaticSimulation sim(1, 1);
        sim.mesh.drawing.drawPolygon(square.copy().translate(center1).rotate(angle1, center1));
        sim.mesh.drawing.drawPolygon(square.copy().translate(center2).rotate(angle2, center2));
        center1.print();
        center2.print();
        sim.mesh.drawing.drawRegion(center1, {0, magnetization1.rotate(angle1), nikfemm::materials::neodymium});
        sim.mesh.drawing.drawRegion(center2, {0, magnetization2.rotate(angle2), nikfemm::materials::neodymium});

        auto system = sim.generateSystem();
        sim.solve(system);

        sim.Bplot(1000, 1000, false, true);

        auto stress1 = sim.computeStressIntegral(center1, center1);
        auto stress2 = sim.computeStressIntegral(center2, center2);

        // F = ma
        // T = Ia
        // update velocity
        velocity1 += stress1.Force * dt / mass1;
        velocity2 += stress2.Force * dt / mass2;
        // update angular velocity
        angular_velocity1 += stress1.Torque * dt / moment_of_inertia1;
        angular_velocity2 += stress2.Torque * dt / moment_of_inertia2;
        // update position
        center1 += velocity1 * dt;
        center2 += velocity2 * dt;
        // update angle
        angle1 += angular_velocity1 * dt;
        angle2 += angular_velocity2 * dt;

        printf("------------------------------------------------------------\n");
        printf("Time: %d\n", i*dt);
        printf("Force 1: %.17g, %.17g\n", stress1.Force.x, stress1.Force.y);
        printf("Force 2: %.17g, %.17g\n", stress2.Force.x, stress2.Force.y);
        printf("Torque 1: %.17g\n", stress1.Torque);
        printf("Torque 2: %.17g\n", stress2.Torque);
        printf("Velocity 1: %.17g, %.17g\n", velocity1.x, velocity1.y);
        printf("Velocity 2: %.17g, %.17g\n", velocity2.x, velocity2.y);
        printf("Angular velocity 1: %.17g\n", angular_velocity1);
        printf("Angular velocity 2: %.17g\n", angular_velocity2);
        printf("Displacement 1: %.17g, %.17g\n", center1.x, center1.y);
        printf("Displacement 2: %.17g, %.17g\n", center2.x, center2.y);
        printf("Angle 1: %.17g\n", angle1);
        printf("Angle 2: %.17g\n", angle2);

        // sim.BplotToFile(400, 400, "magnet_dynamic_simulation_" + std::to_string(i) + ".png");
    }

    // auto stress = simulation.computeStressIntegral({0, 1.5}, {0, 1.5});
    // printf("Force: %.17g, %.17g, Torque: %.17g\n", stress.Force.x, stress.Force.y, stress.Torque);
}