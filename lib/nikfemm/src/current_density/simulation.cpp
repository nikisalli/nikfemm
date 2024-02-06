#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <set>

#ifdef NIKFEMM_USE_OPENCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#endif

#include "../constants.hpp"
#include "simulation.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/vector.hpp"
#include "../algebra/build_coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    CurrentDensitySimulation::CurrentDensitySimulation(double depth) {
        mesh.depth = depth;
    }

    CurrentDensitySimulation::CurrentDensitySimulation() {
        mesh.depth = 1;
        mesh = CurrentDensityMesh();
    }

    CurrentDensitySimulation::~CurrentDensitySimulation() {
        
    }

    void CurrentDensitySimulation::setVoltage(CurrentDensitySystem& system, Vector p, double V) {
        // find the closest node
        int32_t closest_node = -1;
        double closest_distance = INFINITY;
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            double distance = (mesh.data.pointlist[i].x - p.x) * (mesh.data.pointlist[i].x - p.x) + (mesh.data.pointlist[i].y - p.y) * (mesh.data.pointlist[i].y - p.y);
            if (distance < closest_distance) {
                closest_node = i;
                closest_distance = distance;
            }
        }

        if (closest_node == -1) {
            nlogerror("could not find closest node");
            return;
        }

        // set the voltage
        system.addDirichletBoundaryCondition(closest_node, V);
    }

#ifdef NIKFEMM_USE_OPENCV
    void CurrentDensitySimulation::VplotRend(cv::Mat* image, double width, double height) {
        float min_x = mesh.data.pointlist[0].x;
        float min_y = mesh.data.pointlist[0].y;
        float max_x = mesh.data.pointlist[0].x;
        float max_y = mesh.data.pointlist[0].y;
        for (uint32_t i = 1; i < mesh.data.numberofpoints; i++) {
            if (mesh.data.pointlist[i].x < min_x) {
                min_x = mesh.data.pointlist[i].x;
            }
            if (mesh.data.pointlist[i].y < min_y) {
                min_y = mesh.data.pointlist[i].y;
            }
            if (mesh.data.pointlist[i].x > max_x) {
                max_x = mesh.data.pointlist[i].x;
            }
            if (mesh.data.pointlist[i].y > max_y) {
                max_y = mesh.data.pointlist[i].y;
            }
        }

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * width / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * height / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * width - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * height - 0.5 * (max_y + min_y) * y_scale;

        std::vector<float> V_sorted(V.val.size());
        std::copy(V.val.begin(), V.val.end(), V_sorted.begin());
        std::sort(V_sorted.begin(), V_sorted.end());
        float max_V = V_sorted[0.9 * V_sorted.size()];
        float min_V = V_sorted[0.1 * V_sorted.size()];

        nloginfo("max V: %f", max_V);
        nloginfo("min V: %f", min_V);

        // draw the points
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector p = mesh.data.pointlist[i];

            auto points = std::vector<cv::Point>();

            cv::Scalar c = val2jet(V[i], min_V, max_V);

            for (int32_t j = 0; j < mesh.data.numberoftriangles; j++) {
                if (mesh.data.trianglelist[j][0] == i || mesh.data.trianglelist[j][1] == i || mesh.data.trianglelist[j][2] == i) {
                    Vector barycenter = {
                        (mesh.data.pointlist[mesh.data.trianglelist[j][0]].x + mesh.data.pointlist[mesh.data.trianglelist[j][1]].x + mesh.data.pointlist[mesh.data.trianglelist[j][2]].x) / 3,
                        (mesh.data.pointlist[mesh.data.trianglelist[j][0]].y + mesh.data.pointlist[mesh.data.trianglelist[j][1]].y + mesh.data.pointlist[mesh.data.trianglelist[j][2]].y) / 3
                    };
                    points.push_back(cv::Point(x_offset + barycenter.x * x_scale, y_offset + barycenter.y * y_scale));
                }
            }

            // find the center of the points
            cv::Point center = {0, 0};
            for (uint8_t j = 0; j < points.size(); j++) {
                center.x += points[j].x;
                center.y += points[j].y;
            }
            center.x /= points.size();
            center.y /= points.size();

            // sort the points by angle
            std::sort(points.begin(), points.end(), [center](cv::Point a, cv::Point b) {
                return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
            });

            // draw the polygon
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, points), c);
        }
    }

    void CurrentDensitySimulation::Vplot(uint32_t width, uint32_t height) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
        
        if (!image.data) {
            nexit("Could not create image");
        }

        VplotRend(&image, width, height);

        // flip image horizontally
        cv::flip(image, image, 0);

        cv::imshow("V", image);
        nloginfo("showing image");
        cv::waitKey(0);
        nloginfo("done");
    }

    void CurrentDensitySimulation::ElemScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        // get mesh enclosing rectangle
        float min_x = mesh.data.pointlist[0].x;
        float min_y = mesh.data.pointlist[0].y;
        float max_x = mesh.data.pointlist[0].x;
        float max_y = mesh.data.pointlist[0].y;
        for (uint32_t i = 1; i < mesh.data.numberofpoints; i++) {
            if (mesh.data.pointlist[i].x < min_x) {
                min_x = mesh.data.pointlist[i].x;
            }
            if (mesh.data.pointlist[i].y < min_y) {
                min_y = mesh.data.pointlist[i].y;
            }
            if (mesh.data.pointlist[i].x > max_x) {
                max_x = mesh.data.pointlist[i].x;
            }
            if (mesh.data.pointlist[i].y > max_y) {
                max_y = mesh.data.pointlist[i].y;
            }
        }

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * width / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * height / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * width - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * height - 0.5 * (max_y + min_y) * y_scale;

        // std::sort(scalar.begin(), scalar.end());
        // double max_Scalar = scalar[scalar.size() - 1];
        // double min_Scalar = scalar[0];
        std::vector<double> sortedScalar = scalar;
        std::sort(sortedScalar.begin(), sortedScalar.end());
        // 95th percentile
        double max_Scalar = sortedScalar[sortedScalar.size() * 0.95];
        double min_Scalar = sortedScalar[sortedScalar.size() * 0.05];

        nloginfo("max_Scalar: %f, min_Scalar: %f", max_Scalar, min_Scalar);
        

        // max_B = 1e-6;
        // make background grey
        image->setTo(cv::Scalar(127, 127, 127));

        // render
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Vector v1 = mesh.data.pointlist[e[0]];
            Vector v2 = mesh.data.pointlist[e[1]];
            Vector v3 = mesh.data.pointlist[e[2]];

            // get the B vectors
            cv::Scalar c = val2jet(scalar[i], min_Scalar, max_Scalar);
            // get the vertices in window coordinates
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, std::vector<cv::Point>({
                cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset)
            })), c);

            // draw a small white square at each vertex
            cv::rectangle(*image, cv::Point(x_scale * static_cast<float>(v1.x) + x_offset - 1, y_scale * static_cast<float>(v1.y) + y_offset - 1),
                            cv::Point(x_scale * static_cast<float>(v1.x) + x_offset + 1, y_scale * static_cast<float>(v1.y) + y_offset + 1),
                            cv::Scalar(255, 255, 255), 1);
            cv::rectangle(*image, cv::Point(x_scale * static_cast<float>(v2.x) + x_offset - 1, y_scale * static_cast<float>(v2.y) + y_offset - 1),
                            cv::Point(x_scale * static_cast<float>(v2.x) + x_offset + 1, y_scale * static_cast<float>(v2.y) + y_offset + 1),
                            cv::Scalar(255, 255, 255), 1);
            cv::rectangle(*image, cv::Point(x_scale * static_cast<float>(v3.x) + x_offset - 1, y_scale * static_cast<float>(v3.y) + y_offset - 1),
                            cv::Point(x_scale * static_cast<float>(v3.x) + x_offset + 1, y_scale * static_cast<float>(v3.y) + y_offset + 1),
                            cv::Scalar(255, 255, 255), 1);

            // draw mesh edges
            // draw lines from v1 to v2
            if (plotMesh) {
                cv::line(*image, cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                                cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                                cv::Scalar(0, 0, 0), 1);
                cv::line(*image, cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                                cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset),
                                cv::Scalar(0, 0, 0), 1);
                cv::line(*image, cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset),
                                cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                                cv::Scalar(0, 0, 0), 1);
            }
        }

        // draw the geometry
        // draw the segments
        for (DrawingSegment s : mesh.drawing.segments) {
            cv::line(*image, cv::Point(x_scale * mesh.drawing.points[s.p1].x + x_offset,
                                       y_scale * mesh.drawing.points[s.p1].y + y_offset), 
                             cv::Point(x_scale * mesh.drawing.points[s.p2].x + x_offset, 
                                       y_scale * mesh.drawing.points[s.p2].y + y_offset),
                             cv::Scalar(255, 255, 255), 1);
        }

        if (plotRegions) {
            // draw the regions
            for (DrawingRegion r : mesh.drawing.regions) {
                std::vector<cv::Point> points;
                
                Vector pos = r.first;
                // draw white cross
                cv::line(*image, cv::Point(x_scale * pos.x + x_offset - 5, y_scale * pos.y + y_offset),
                                 cv::Point(x_scale * pos.x + x_offset + 5, y_scale * pos.y + y_offset),
                                 cv::Scalar(255, 255, 255), 1);
                cv::line(*image, cv::Point(x_scale * pos.x + x_offset, y_scale * pos.y + y_offset - 5),
                                 cv::Point(x_scale * pos.x + x_offset, y_scale * pos.y + y_offset + 5),
                                 cv::Scalar(255, 255, 255), 1);
            }
        }
    }

    void CurrentDensitySimulation::ElemScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        ElemScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        // flip image horizontally
        cv::flip(image, image, 0);

        // show the image
        cv::imshow("B", image);
        // continue if image is closed
        cv::waitKey(0);
    }
#endif

    CurrentDensitySystem CurrentDensitySimulation::generateSystem(bool refine, double max_triangle_area, int min_angle) {
        if (refine) {
            mesh.drawing.addRefiningPoints();
        }

        mesh.mesh(max_triangle_area, min_angle);
        mesh.computeEpsilon();
        nloginfo("the mesh has %u nodes and %u elements", mesh.data.numberofpoints, mesh.data.numberoftriangles);

        auto system = mesh.getFemSystem();
        // auto system = mesh.getFemSystemCotangentWeights();

        return system;
    }

    void CurrentDensitySimulation::solve(CurrentDensitySystem& system) {
        V = CV(mesh.data.numberofpoints);

        MatCSRSymmetric FemMat(system.A);

        auto start = std::chrono::high_resolution_clock::now();
        preconditionedSSORConjugateGradientSolver(FemMat, system.b, V, 1.5, 1e-12, 100000);
        // preconditionedJacobiConjugateGradientSolver(FemMat, system.b, V, 1e-6, 100000);
        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    }
}