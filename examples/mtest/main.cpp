#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <atomic>
#include <cassert>
#include <chrono>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

using namespace nikfemm;

// too lazy to write my own thread pool
// https://maidamai0.github.io/post/a-simple-thread-pool/
class ThreadPool {
    using task_type = std:: function < void() > ;

    public:
        ThreadPool(size_t num = std::thread::hardware_concurrency()) {
            for (size_t i = 0; i < num; ++i) {
                workers_.emplace_back(std::thread([this] {
                    while (true) {
                        task_type task; {
                            std::unique_lock < std::mutex > lock(task_mutex_);
                            task_cond_.wait(lock, [this] {
                                return !tasks_.empty();
                            });
                            task = std::move(tasks_.front());
                            tasks_.pop();
                        }
                        if (!task) {
                            std::cout << "worker #" << std::this_thread::get_id() << " exited" << std::endl;
                            push_stop_task();
                            return;
                        }
                        task();
                    }
                }));
                std::cout << "worker #" << workers_.back().get_id() << " started" << std::endl;
            }
        }

        ~ThreadPool() {
            Stop();
        }

    void Stop() {
        push_stop_task();
        for (auto & worker: workers_) {
            if (worker.joinable()) {
                worker.join();
            }
        }

        // clear all pending tasks
        std::queue < task_type > empty {};
        std::swap(tasks_, empty);
    }

    template < typename F, typename...Args >
        auto Push(F && f, Args && ...args) {
            using return_type = std::invoke_result_t < F, Args... > ;
            auto task
                = std::make_shared < std::packaged_task < return_type() >> (std::bind(std::forward < F > (f), std::forward < Args > (args)...));
            auto res = task -> get_future();

            {
                std::lock_guard < std::mutex > lock(task_mutex_);
                tasks_.emplace([task] {
                    ( * task)();
                });
            }
            task_cond_.notify_one();

            return res;
        }

    private:
        void push_stop_task() {
            std::lock_guard < std::mutex > lock(task_mutex_);
            tasks_.push(task_type {});
            task_cond_.notify_one();
        }

    std::vector < std::thread > workers_;
    std::queue < task_type > tasks_;
    std::mutex task_mutex_;
    std::condition_variable task_cond_;
};

Vector polar_to_scalar(double r, double theta) {
    return Vector(r * cos(theta), r * sin(theta));
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

    std::vector<Vector> stator;
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
        // printf("winding %d inverted %d\n", winding_id, i % 2);

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

    sim.mesh.drawing.drawCircle(Vector(0, 0), stator_inner_radius, stator_poles * 10);

    // stator core
    sim.mesh.drawing.drawRegion(polar_to_scalar(stator_outer_radius, 0), {0, {0, 0}, nikfemm::materials::iron_linear});

    // inter stator-rotor gap
    sim.mesh.drawing.drawRegion(polar_to_scalar(
        stator_outer_radius + winding_to_pole_and_core_distance / 2,
        angle_step / 2
    ), {0, {0, 0}, nikfemm::materials::air});

    // stator center hole
    sim.mesh.drawing.drawRegion(Vector(0, 0), {0, {0, 0}, nikfemm::materials::air});

    // draw rotor
    sim.mesh.drawing.drawCircle(Vector(0, 0), rotor_outer_radius, rotor_poles * 10);

    // draw magnets
    double magnet_angle_step = (2 * PI) / (rotor_poles);
    double magnet_width_half = magnet_angle_step * magnet_fill_width_percent / 2;

    std::vector<Vector> magnet_inner_rotor;

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
            {0, polar_to_scalar(magnet_orientation, center_angle), nikfemm::materials::iron_linear}
        );
        
        magnet_inner_rotor.push_back(polar_to_scalar(rotor_inner_radius, center_angle - magnet_width_half));
        magnet_inner_rotor.push_back(polar_to_scalar(rotor_inner_radius, center_angle + magnet_width_half));
    }

    sim.mesh.drawing.drawPolygon(magnet_inner_rotor);
    // rotor
    sim.mesh.drawing.drawRegion(polar_to_scalar((rotor_outer_radius + rotor_inner_radius) / 2, angle), {0, {0, 0}, nikfemm::materials::iron_linear});
}

void simulate(int num, double angle, double max_B, double min_B,
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
    MagnetostaticSimulation simulation;

    draw_motor(simulation,
        angle,
        stator_poles,
        rotor_poles,
        stator_outer_pole_radius,
        stator_inner_pole_radius,
        stator_outer_radius,
        stator_inner_radius,
        stator_pole_cap_percent_unfilled,
        stator_pole_stem_percent_of_angle_step,
        rotor_outer_radius,
        rotor_inner_radius,
        magnet_thickness,
        magnet_fill_width_percent,
        winding_fill_width_percent,
        winding_percent_of_stem,
        currentA,
        currentB,
        currentC,
        magnet_magnetization
    );

    simulation.solve();
    // simulation.AplotToFile(10000, 10000, "Aplot.png");
    // save to file num
    simulation.BplotToFile(10000, 10000, "Bplot" + std::to_string(num) + ".png", false, true);
}

int main(int argc, char** argv) {
    int num = 1;
    double angle = 0;
    // calculate max and min B by simulating at 0 angle
    double max_B = 0.000073;
    double min_B = 0.000000;
    // MagnetostaticSimulation simulation;
    // draw_motor(simulation, 0, 24, 20, 1, 0.95, 0.45, 0.30, 0.1, 0.5, 1.13, 1.12, 0.10, 0.9, 0.24, 0.95, 0, 0.8660254038, -0.8660254038, 1);
    // simulation.solve();
    // std::vector<double> B_mags;
    // for (uint32_t i = 0; i < simulation.B.size(); i++) {
    //     B_mags.push_back(simulation.B[i].magnitude());
    // }
    // std::sort(B_mags.begin(), B_mags.end());
    // max_B = B_mags[0.9 * B_mags.size()];
    // min_B = B_mags[0.1 * B_mags.size()];
    // printf("max B: %f, min B: %f\n", max_B, min_B);
    
    ThreadPool pool(12);
    for (double angle = 0; angle < (2 * PI) / 4; angle += ((2 * PI) / 4) / 1000) {
        pool.Push([angle, num, max_B, min_B]() {
            simulate(num, angle, max_B, min_B,
                24, // stator_poles
                20, // rotor_poles
                1, // stator_outer_pole_radius
                0.95, // stator_inner_pole_radius
                0.45, // stator_outer_radius
                0.30, // stator_inner_radius
                0.1, // stator_pole_cap_percent_unfilled
                0.5, // stator_pole_stem_percent_of_angle_step
                1.20, // rotor_outer_radius
                1.12, // rotor_inner_radius
                0.10, // magnet_thickness
                0.7, // magnet_fill_width_percent
                0.24, // winding_fill_width_percent
                0.95, // winding_percent_of_stem
                sin(angle), // currentA
                sin(angle + (2 * PI) / 3), // currentB
                sin(angle + (4 * PI) / 3), // currentC
                0.5 // magnet_magnetization
            );
        });
        num++;
    }
    return 0;
}