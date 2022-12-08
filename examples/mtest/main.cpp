#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace nikfemm;

Point polar_to_scalar(double r, double theta) {
    return Point(r * cos(theta), r * sin(theta));
}

void draw_motor(MagnetostaticSimulation& sim, double angle, 
                int stator_poles,
                int rotor_poles,
                double stator_outer_pole_radius,
                double stator_inner_pole_radius,
                double stator_outer_radius, 
                double stator_inner_radius,
                double stator_pole_cap_percent_unfilled,
                double stator_pole_stem_percent_of_angle_step,
                double rotor_outer_radius,
                double rotor_inner_radius,
                double magnet_thickness,
                double magnet_fill_width_percent,
                double winding_fill_width_percent,
                double winding_percent_of_stem,
                double currentA,
                double currentB,
                double currentC,
                double magnet_magnetization) {

    // draw stator
    double angle_step = (2 * PI) / (stator_poles);
    double pole_stem_half = angle_step * stator_pole_stem_percent_of_angle_step / 2;
    double pole_cap_half = angle_step * (1 - stator_pole_cap_percent_unfilled) / 2;
    double winding_angle_step = angle_step * winding_fill_width_percent;
    double stem_height = stator_inner_pole_radius - stator_outer_radius;
    double winding_height = stem_height * winding_percent_of_stem;
    double winding_to_pole_and_core_distance = stem_height - winding_height;
    double winding_start_radius = stator_outer_radius + winding_to_pole_and_core_distance;
    double winding_end_radius = stator_inner_pole_radius - winding_to_pole_and_core_distance;

    std::vector<Point> stator;
    for (int i = 0; i < stator_poles; i++) {
        double center_angle = angle_step * i;
        stator.push_back(polar_to_scalar( stator_outer_radius, center_angle - pole_stem_half));
        stator.push_back(polar_to_scalar( winding_start_radius, center_angle - pole_stem_half));
        stator.push_back(polar_to_scalar( winding_end_radius, center_angle - pole_stem_half));
        stator.push_back(polar_to_scalar( stator_inner_pole_radius, center_angle - pole_stem_half));
        stator.push_back(polar_to_scalar( stator_inner_pole_radius, center_angle - pole_cap_half));
        stator.push_back(polar_to_scalar( stator_outer_pole_radius, center_angle - pole_cap_half));
        stator.push_back(polar_to_scalar( stator_outer_pole_radius, center_angle + pole_cap_half));
        stator.push_back(polar_to_scalar( stator_inner_pole_radius, center_angle + pole_cap_half));
        stator.push_back(polar_to_scalar( stator_inner_pole_radius, center_angle + pole_stem_half));
        stator.push_back(polar_to_scalar( winding_end_radius, center_angle + pole_stem_half));
        stator.push_back(polar_to_scalar( winding_start_radius, center_angle + pole_stem_half));
        stator.push_back(polar_to_scalar( stator_outer_radius, center_angle + pole_stem_half));

        // draw winding
        int winding_id = (i / 2) % 3;
        float current = 0;
        switch (winding_id) {
            case 0:
                current = currentA;
                break;
            case 1:
                current = currentB;
                break;
            case 2:
                current = currentC;
                break;
        }

        if ((i / 2) % 2) current *= -1;
        printf("winding %d inverted %d\n", winding_id, i % 2);

        sim.mesh.drawing.drawPolygon({
            polar_to_scalar(winding_start_radius, center_angle - pole_stem_half),
            polar_to_scalar(winding_start_radius, center_angle - pole_stem_half - winding_angle_step),
            polar_to_scalar(winding_end_radius, center_angle - pole_stem_half - winding_angle_step),
            polar_to_scalar(winding_end_radius, center_angle - pole_stem_half),
        });
        // winding
        sim.mesh.drawing.drawRegion(polar_to_scalar(
            (winding_end_radius - winding_start_radius) / 2 + winding_start_radius,
            center_angle - pole_stem_half - winding_angle_step / 2
        ), {-current, {0, 0}, nikfemm::materials::air});

        sim.mesh.drawing.drawPolygon({
            polar_to_scalar(winding_start_radius, center_angle + pole_stem_half),
            polar_to_scalar(winding_start_radius, center_angle + pole_stem_half + winding_angle_step),
            polar_to_scalar(winding_end_radius, center_angle + pole_stem_half + winding_angle_step),
            polar_to_scalar(winding_end_radius, center_angle + pole_stem_half),
        });
        // winding
        sim.mesh.drawing.drawRegion(polar_to_scalar(
            (winding_end_radius - winding_start_radius) / 2 + winding_start_radius,
            center_angle + pole_stem_half + winding_angle_step / 2
        ), {current, {0, 0}, nikfemm::materials::air});

        // draw air region
    }

    sim.mesh.drawing.drawPolygon(stator);

    sim.mesh.drawing.drawCircle(Point(0, 0), stator_inner_radius, stator_poles * 10);

    // stator core
    sim.mesh.drawing.drawRegion(polar_to_scalar(stator_outer_radius, 0), {0, {0, 0}, nikfemm::materials::iron_linear});

    // inter stator-rotor gap
    sim.mesh.drawing.drawRegion(polar_to_scalar(
        stator_outer_radius + winding_to_pole_and_core_distance / 2,
        angle_step / 2
    ), {0, {0, 0}, nikfemm::materials::air});

    // stator center hole
    sim.mesh.drawing.drawRegion(Point(0, 0), {0, {0, 0}, nikfemm::materials::air});

    // draw rotor
    sim.mesh.drawing.drawCircle(Point(0, 0), rotor_outer_radius, rotor_poles * 10);

    // draw magnets
    double magnet_angle_step = (2 * PI) / (rotor_poles);
    double magnet_width_half = magnet_angle_step * magnet_fill_width_percent / 2;

    std::vector<Point> magnet_inner_rotor;

    for (int i = 0; i < rotor_poles; i++) {
        double center_angle = magnet_angle_step * i + angle;

        sim.mesh.drawing.drawPolygon({
            polar_to_scalar(rotor_inner_radius, center_angle - magnet_width_half),
            polar_to_scalar(rotor_inner_radius - magnet_thickness, center_angle - magnet_width_half),
            polar_to_scalar(rotor_inner_radius - magnet_thickness, center_angle + magnet_width_half),
            polar_to_scalar(rotor_inner_radius, center_angle + magnet_width_half)
        });

        double magnet_orientation = (i % 2 == 0) ? -magnet_magnetization : magnet_magnetization;
        // magnet
        sim.mesh.drawing.drawRegion(
            polar_to_scalar(rotor_inner_radius - (magnet_thickness / 2), center_angle),
            {0, polar_to_scalar(magnet_orientation, center_angle), nikfemm::materials::neodymium}
        );
        
        magnet_inner_rotor.push_back(polar_to_scalar(rotor_inner_radius, center_angle - magnet_width_half));
        magnet_inner_rotor.push_back(polar_to_scalar(rotor_inner_radius, center_angle + magnet_width_half));
    }

    sim.mesh.drawing.drawPolygon(magnet_inner_rotor);
    // rotor
    sim.mesh.drawing.drawRegion(polar_to_scalar((rotor_outer_radius + rotor_inner_radius) / 2, angle), {0, {0, 0}, nikfemm::materials::iron_linear});
}

int main(int argc, char** argv) {
    int num = 1;
    double angle = 0;
    // for (double angle = 0; angle < (2 * PI) / 4; angle += ((2 * PI) / 4) / 1000) {
    // printf("angle %f\n", angle);
    MagnetostaticSimulation simulation;

    draw_motor(simulation,
                angle, // angle
                24, // stator_poles
                20, // rotor_poles
                1, // stator_outer_pole_radius
                0.95, // stator_inner_pole_radius
                0.45, // stator_outer_radius
                0.30, // stator_inner_radius
                0.1, // stator_pole_cap_percent_unfilled
                0.5, // stator_pole_stem_percent_of_angle_step
                1.40, // rotor_outer_radius
                1.12, // rotor_inner_radius
                0.10, // magnet_thickness
                0.9, // magnet_fill_width_percent
                0.24, // winding_fill_width_percent
                0.95, // winding_percent_of_stem
                1.0, // currentA
                0.866025, // currentB
                -0.866025, // currentC
                1 // magnet_magnetization

    );

    simulation.solve();
    // simulation.AplotToFile(10000, 10000, "Aplot.png");
    // save to file num
    simulation.BplotToFile(10000, 10000, "Bplot" + std::to_string(num) + ".png", true, false);
        // num++;
    // }
    return 0;
}