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

#include <omp.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>

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
    MagnetostaticSimulation::MagnetostaticSimulation() {

    }

    MagnetostaticSimulation::~MagnetostaticSimulation() {

    }

    void MagnetostaticSimulation::AplotRend(cv::Mat* image, double width, double height) {
        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

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

        // find max and min
        /*
        float max_A = std::numeric_limits<float>::min();
        float min_A = std::numeric_limits<float>::max();

        for (uint32_t i = 0; i < A.m; i++) {
            if (A[i] > max_A) {
                max_A = A[i];
            }
            if (A[i] < min_A) {
                min_A = A[i];
            }
        }
        */

        std::vector<float> A_sorted(A.val.size());
        for (uint32_t i = 0; i < A.val.size(); i++) {
            A_sorted[i] = A[i];
        }
        std::sort(A_sorted.begin(), A_sorted.end());
        float max_A = A_sorted[0.9 * A_sorted.size()];
        float min_A = A_sorted[0.1 * A_sorted.size()];

        printf("max_A: %f, min_A: %f\n", max_A, min_A);

        // draw the points
        for (int64_t i = 0; i < mesh.data.numberofpoints; i++) {    
            Vector p = mesh.data.pointlist[i];

            if (Vector::distance(p, Vector(0, 0)) > mesh.radius + mesh.epsilon) {
                continue;
            }

            // verts
            auto points = std::vector<cv::Point>();

            // fill the Vertex<MagnetostaticProp> cell with jet color of the Vertex<MagnetostaticProp>
            cv::Scalar c = val2jet(A[i], min_A, max_A);
            // printf("A: %f, c: %d, %d, %d\n", A.val[i], c.r, c.g, c.b);
                            
            // find the triangles that contain the Vertex<MagnetostaticProp> and then
            // for every triangle find the barycenter and add it to the points vector
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
            // printf("count: %d\n", points.size());
            center.x /= points.size();
            center.y /= points.size();
            // printf("center: %d %d\n", center.x, center.y);

            // sort the points by angle
            std::sort(points.begin(), points.end(), [center](cv::Point a, cv::Point b) {
                return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
            });

            // draw the polygon
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, points), c);
        }
    }

    void MagnetostaticSimulation::Aplot(uint32_t width, uint32_t height) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        if (!image.data) {
            nexit("Could not create image");
        }

        AplotRend(&image, width, height);

        cv::imshow("A", image);
        printf("showing image\n");
        cv::waitKey(0);
        printf("done\n");
    }

    void MagnetostaticSimulation::AplotToFile(uint32_t width, uint32_t height, std::string filename) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        if (!image.data) {
            nexit("Could not create image");
        }

        AplotRend(&image, width, height);

        cv::imwrite(filename, image);
    }

    void MagnetostaticSimulation::NodeScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

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

        std::vector<float> A_sorted(scalar.size());
        for (uint32_t i = 0; i < scalar.size(); i++) {
            A_sorted[i] = scalar[i];
        }
        std::sort(A_sorted.begin(), A_sorted.end());
        float max_A = A_sorted[0.9 * A_sorted.size()];
        float min_A = A_sorted[0.1 * A_sorted.size()];

        printf("max_A: %f, min_A: %f\n", max_A, min_A);

        // draw the points
        for (int64_t i = 0; i < mesh.data.numberofpoints; i++) {    
            Vector p = mesh.data.pointlist[i];

            if (Vector::distance(p, Vector(0, 0)) > mesh.radius * 1.5 + mesh.epsilon) {
                continue;
            }

            // verts
            auto points = std::vector<cv::Point>();

            // fill the Vertex<MagnetostaticProp> cell with jet color of the Vertex<MagnetostaticProp>
            cv::Scalar c = val2jet(scalar[i], min_A, max_A);
            // printf("A: %f, c: %d, %d, %d\n", A.val[i], c.r, c.g, c.b);
                            
            // find the triangles that contain the Vertex<MagnetostaticProp> and then
            // for every triangle find the barycenter and add it to the points vector
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
            // printf("count: %d\n", points.size());
            center.x /= points.size();
            center.y /= points.size();
            // printf("center: %d %d\n", center.x, center.y);

            // sort the points by angle
            std::sort(points.begin(), points.end(), [center](cv::Point a, cv::Point b) {
                return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
            });

            // draw the polygon
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, points), c);
        }
    }

    void MagnetostaticSimulation::NodeScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        if (!image.data) {
            nexit("Could not create image");
        }

        NodeScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        cv::imshow("A", image);
        printf("showing image\n");
        cv::waitKey(0);
        printf("done\n");
    }

    void MagnetostaticSimulation::NodeScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh, bool plotRegions) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        if (!image.data) {
            nexit("Could not create image");
        }

        NodeScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        cv::imwrite(filename, image);
    }

    void MagnetostaticSimulation::BplotRend(cv::Mat* image, double width, double height, bool plotMesh, bool plotRegions, double maxB, double minB, bool plotCurves, uint32_t curve_number) {
        // get mesh enclosing rectangle
        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

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

        // find max and min B
        /*
        double max_B = std::numeric_limits<float>::min();
        double min_B = std::numeric_limits<float>::max();

        for (uint32_t i = 0; i < B.size(); i++) {
            double B_mag = B[i].magnitude();
            if (B_mag > max_B) {
                max_B = B_mag;
            }
            if (B_mag < min_B) {
                min_B = B_mag;
            }
        }
        */

        double max_B;
        double min_B;

        if (std::isnan(maxB) || std::isnan(minB)) {
            std::vector<double> B_mags;
            for (uint32_t i = 0; i < B.size(); i++) {
                B_mags.push_back(B[i].norm());
            }
            std::sort(B_mags.begin(), B_mags.end());
            max_B = B_mags[0.99 * B_mags.size()];
            min_B = B_mags[0.01 * B_mags.size()];
        } else {
            max_B = maxB;
            min_B = minB;
        }

        printf("max B: %f\n", max_B);
        printf("min B: %f\n", min_B);

        // max_B = 1e-6;

        // render
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Vector v1 = mesh.data.pointlist[e[0]];
            Vector v2 = mesh.data.pointlist[e[1]];
            Vector v3 = mesh.data.pointlist[e[2]];

            if (Vector::distance(v1, Vector(0, 0)) > mesh.radius + mesh.epsilon ||
                Vector::distance(v2, Vector(0, 0)) > mesh.radius + mesh.epsilon || 
                Vector::distance(v3, Vector(0, 0)) > mesh.radius + mesh.epsilon) {
                continue;
            }

            // get the B vectors
            Vector Bv = B[i];
            cv::Scalar c = val2jet(Bv.norm(), min_B, max_B);
            // get the vertices in window coordinates
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, std::vector<cv::Point>({
                cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset)
            })), c);

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

            // level curves
            
        }

        // draw the geometry
        // draw the segments
        for (DrawingSegment s : mesh.drawing.segments) {
            cv::line(*image, cv::Point(x_scale * mesh.drawing.points[s.p1].x + x_offset,
                                       y_scale * mesh.drawing.points[s.p1].y + y_offset), 
                             cv::Point(x_scale * mesh.drawing.points[s.p2].x + x_offset, 
                                       y_scale * mesh.drawing.points[s.p2].y + y_offset),
                             cv::Scalar(0, 0, 0), 1);
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
                // draw an arrow for each magnet region
                /*
                */
                if (mesh.drawing.region_map[r.second].M != Vector(0, 0)) {
                    Vector start = r.first;
                    Vector end = start + mesh.drawing.region_map[r.second].M;
                    cv::arrowedLine(*image, cv::Point(x_scale * start.x + x_offset, y_scale * start.y + y_offset),
                                            cv::Point(x_scale * end.x + x_offset, y_scale * end.y + y_offset),
                                            cv::Scalar(255, 255, 255), 1);
                }
            }

            // draw all drawing points
            for (Vector p : mesh.drawing.points) {
                cv::circle(*image, cv::Point(x_scale * p.x + x_offset, y_scale * p.y + y_offset), 2, cv::Scalar(255, 255, 255), -1);
            }
        }
    }

    void MagnetostaticSimulation::Bplot(uint32_t width, uint32_t height, bool plotMesh, bool plotRegions, double maxB, double minB, bool plotCurves, uint32_t curve_number) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        BplotRend(&image, width, height, plotMesh, plotRegions, maxB, minB, plotCurves, curve_number);

        // show the image
        cv::imshow("B", image);
        // continue if image is closed
        cv::waitKey(0);
    }

    void MagnetostaticSimulation::BplotToFile(uint32_t width, uint32_t height, std::string filename, bool plotMesh, bool plotRegions, double maxB, double minB, bool plotCurves, uint32_t curve_number) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        BplotRend(&image, width, height, plotMesh, plotRegions, maxB, minB, plotCurves, curve_number);

        // save the image
        cv::imwrite(filename, image);
    }

    void MagnetostaticSimulation::ElemScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        // get mesh enclosing rectangle
        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

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

        printf("max_Scalar: %f, min_Scalar: %f\n", max_Scalar, min_Scalar);

        // max_B = 1e-6;

        // render
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Vector v1 = mesh.data.pointlist[e[0]];
            Vector v2 = mesh.data.pointlist[e[1]];
            Vector v3 = mesh.data.pointlist[e[2]];

            if (Vector::distance(v1, Vector(0, 0)) > mesh.radius * 1.5 + mesh.epsilon ||
                Vector::distance(v2, Vector(0, 0)) > mesh.radius * 1.5 + mesh.epsilon || 
                Vector::distance(v3, Vector(0, 0)) > mesh.radius * 1.5 + mesh.epsilon) {
                continue;
            }

            // get the B vectors
            cv::Scalar c = val2jet(scalar[i], min_Scalar, max_Scalar);
            // get the vertices in window coordinates
            cv::fillPoly(*image, std::vector<std::vector<cv::Point>>(1, std::vector<cv::Point>({
                cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset)
            })), c);

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

    void MagnetostaticSimulation::ElemScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        ElemScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        // show the image
        cv::imshow("B", image);
        // continue if image is closed
        cv::waitKey(0);
    }

    void MagnetostaticSimulation::ElemScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        ElemScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        // save the image
        cv::imwrite(filename, image);
    }

    void MagnetostaticSimulation::updateMu(std::vector<const MagnetostaticProp*>& props, std::vector<float>& mu, std::vector<Vector>& B, double residual, uint32_t iter) {
        assert(mu.size() == B.size());
        double kalman_scemo = exp(-((double)iter) / 5.);
        // printf("kalman_scemo: %.17g\n", kalman_scemo);
        // printf("exp(-(iter + 1) / 100): %.17g\n", exp(-((double)iter + 1.) / 100.));
        // printf("- (iter + 1) / 100: %.17g\n", - ((double)iter + 1.) / 100.);
        
        for (uint32_t i = 0; i < B.size(); i++) {
            float Bmag = B[i].norm();
            if (props[i]->isLinear()) {
                mu[i] = props[i]->getMu(Bmag);
            } else {
                mu[i] = (props[i]->getMu(Bmag) * kalman_scemo) + ((mu[i] + materials::vacuum * residual * residual) * (1 - kalman_scemo)); 
                // mu[i] += (props[i]->getMu(Bmag) - mu[i]) * 0.1;  // mu += (mu_new - mu) * 0.1
                // mu[i] += props[i]->getMu(Bmag);  // pure newton-raphson
            }
        }
    }

    MagnetostaticSystem MagnetostaticSimulation::generateSystem(bool refine) {
        // get time in milliseconds
        if (refine) {
            mesh.drawing.addRefiningPoints();
        }

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(mesh.drawing.points);
        if (mesh.drawing.points.size() == 0) {
            smallest_circle.radius = 1;
        }
        // set simulation offset and boundary radius
        mesh.center = smallest_circle.center;
        mesh.radius = smallest_circle.radius * 2;
        // translate everything to the origin
        mesh.drawing.translate(-mesh.center);
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(Vector(0, 0), 2 * smallest_circle.radius);
        double circumferential_length = boundary_circle.circumference();
        uint32_t boundary_points = (uint32_t)((circumferential_length / sqrt((MAX_TRIANGLE_AREA * 4) / sqrt(3))) * 2);
        printf("boundary points: %d\n", boundary_points);
        mesh.drawing.drawCircle(boundary_circle, boundary_points);
        // add region near the edge of the circle
        mesh.drawing.drawRegion(Vector(boundary_circle.radius * 0.9, 0), {0, {0, 0}, materials::air});
        // add the boundary 
        // mesh.drawing.plot();
        #ifdef NIK_REFINE_MAGNETS
            mesh.refineMeshAroundMagnets();
        #endif
        mesh.mesh();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        mesh.addKelvinBoundaryConditions(boundary_points);
        mesh.computeEpsilon();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        printf("the mesh has %u nodes and %u elements\n", mesh.data.numberofpoints, mesh.data.numberoftriangles);

        auto system = mesh.getFemSystem();
        mesh.addDirichletInfiniteBoundaryConditions(system);
        return system;
    }

    void MagnetostaticSimulation::solve(MagnetostaticSystem& system) {
        A = CV(mesh.data.numberofpoints);
        B = std::vector<Vector>(mesh.data.numberoftriangles, {0, 0});
        std::vector<float> mu(mesh.data.numberoftriangles, 0);
        std::vector<const MagnetostaticProp*> props(mesh.data.numberoftriangles);

        // fill props
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            props[i] = mesh.drawing.getRegionPtrFromId(mesh.data.triangleattributelist[i]);
        }

        MagnetostaticMatCSRSymmetric FemMat(system.coo);

        // initialize mu
        for (uint32_t i = 0; i < B.size(); i++) {
            mu[i] = props[i]->getMu(0);
        }
        FemMat.updateFromMu(mu);
        system.b.updateFromMu(mu);

        // check if materials are all linear
        bool all_linear = true;
        for (uint32_t i = 0; i < props.size(); i++) {
            if (!props[i]->isLinear()) {
                all_linear = false;
                break;
            }
        }

        if (all_linear) {
            printf("all materials are linear\n");
            preconditionedSSORConjugateGradientSolver(FemMat, system.b, A, 1.5, 1e-6, 100000);
            // preconditionedIncompleteCholeskyConjugateGradientSolver(FemMat, b, A, 1e-6, 100000);
        } else {
            printf("nonlinear materials detected, starting non linear newton solver\n");
            CV r(system.b.val.size());  // residual

            double residual = 1e10;
            for (uint32_t i = 0; i < 500; i++) {
                // conjugateGradientSolver(FemMat, b, A, 1e-7, 10000);
                // preconditionedJacobiConjugateGradientSolver(FemMat, b, A, 1e-7, 1000);
                preconditionedSSORConjugateGradientSolver(FemMat, system.b, A, 1.5, 1e-7, 50);
                // preconditionedIncompleteCholeskyConjugateGradientSolver(FemMat, b, A, 1e-7, 1000);
                // mesh.Aplot(A);
                // mesh.Bplot(B);

                mesh.computeCurl(B, A);
                // mesh.Bplot(B);

                // check if the solution is correct
                MagnetostaticSimulation::updateMu(props, mu, B, residual, i);
                FemMat.updateFromMu(mu);
                system.b.updateFromMu(mu);
                CV::mult(r, FemMat, A);
                CV::sub(r, system.b, r);
                residual = CV::norm(r);
                if (residual < 1e-7) {
                    printf("Converged in %d iterations\n", i);
                    break;
                }
                // printf("%.17g, %.17g; ", K, residual);
                printf("nonlinear iteration %d, residual: %.17g\n", i, residual);
                // print mu
                // printf("%f,", residual);
                fflush(stdout);
            }
        }

        mesh.computeCurl(B, A);
    }

    Vector MagnetostaticSimulation::computeForceIntegrals(Vector p) {
        // find the polygon that contains p
        Polygon integration_region;
        Polygon boundary_region;
        std::vector<Polygon> boundary_region_holes;
        Vector translated_p = p - mesh.center;
        std::vector<Polygon> polygons_that_contain_p;
        for (auto polygon : mesh.drawing.polygons) {
            if (polygon.contains(translated_p)) {
                polygons_that_contain_p.push_back(polygon);
            }
        }

        // if there is more than one polygon that contains p then we have to find the polygon that isn't contained by any other polygon that contains p
        if (polygons_that_contain_p.size() == 1) {
            integration_region = polygons_that_contain_p[0];
        } else {
            for (auto polygon : polygons_that_contain_p) {
                bool is_integration_region = true;
                for (auto other_polygon : polygons_that_contain_p) {
                    if (polygon != other_polygon && polygon.contains(other_polygon)) {
                        is_integration_region = false;
                        break;
                    }
                }
                if (is_integration_region) {
                    integration_region = polygon;
                    break;
                }
            }
        }

        // find the region that immediately surrounds the integration region
        std::vector<Polygon> polygons_that_contain_integration_region;
        for (auto polygon : mesh.drawing.polygons) {
            if (polygon.contains(integration_region) && polygon != integration_region) {
                polygons_that_contain_integration_region.push_back(polygon);
            }
        }
        // if there is more than one polygon that contains the integration region then we have to keep
        // removing polygons that aren't contained by any other polygon that contains the integration region
        // until we are left with one polygon
        int count = 0;
        while (polygons_that_contain_integration_region.size() > 1) {
            for (auto polygon : polygons_that_contain_integration_region) {
                bool is_outermost_polygon = true;
                for (auto other_polygon : polygons_that_contain_integration_region) {
                    if (polygon != other_polygon && polygon.contains(other_polygon)) {
                        is_outermost_polygon = false;
                        break;
                    }
                }
                if (is_outermost_polygon) {
                    // remove polygon from polygons_that_contain_integration_region
                    polygons_that_contain_integration_region.erase(
                        std::remove(polygons_that_contain_integration_region.begin(),
                                    polygons_that_contain_integration_region.end(),
                                    polygon),
                        polygons_that_contain_integration_region.end());
                }
            }
            count++;
            if (count > 100) {
                nexit("error: could not find outermost polygon");
            }
        }

        // the boundary region is the polygon that is left
        boundary_region = *polygons_that_contain_integration_region.begin();

        // find the holes in the boundary region (every polygon that is contained by the boundary region but isn't contained by the integration region)
        for (auto polygon : mesh.drawing.polygons) {
            if (boundary_region.contains(polygon) && !integration_region.contains(polygon) && polygon != boundary_region && polygon != integration_region) {
                boundary_region_holes.push_back(polygon);
            }
        }

        auto prop = mesh.drawing.getPolygonProp(integration_region);

        // EGGSHELL WITH LAPLACE PROBLEM

        auto coo = BuildMatCOO<double>(mesh.data.numberofpoints);
        auto b = CV(mesh.data.numberofpoints);
        auto b_dirichlet_mask = std::vector<bool>(mesh.data.numberofpoints, false);
        // since the stiffness matrix is symmetric, this function only computes the upper triangular part

        // COMPUTE FEM WEIGHTS
        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(mesh.data.numberofpoints);
        auto adjelems_count = std::vector<uint8_t>(mesh.data.numberofpoints, 0);

        printf("computing adjelems_ids\n");
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = mesh.data.trianglelist[i][j];
                adjelems_ids[myid][adjelems_count[myid]++] = i;
            }
        }

        printf("computing elemadjelems_ids\n");
        struct Edge {
            uint8_t count = 0;
            uint32_t elems[2];
        };
        uint64_t edge_count = (mesh.data.numberofpoints - 2) * 3;
        std::unordered_map<uint64_t, Edge> edge_to_triangles(edge_count);
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t v1 = mesh.data.trianglelist[i][j];
                uint32_t v2 = mesh.data.trianglelist[i][(j + 1) % 3];
                uint64_t edge = v1 < v2 ? (uint64_t)v1 << 32 | v2 : (uint64_t)v2 << 32 | v1;
                edge_to_triangles[edge].elems[edge_to_triangles[edge].count++] = i;
            }
        }

        printf("computing field error\n");
        // compute field error for each element
        std::vector<double> elem_field_errors = std::vector<double>(mesh.data.numberoftriangles, 0);
        std::vector<double> elem_field_weights = std::vector<double>(mesh.data.numberoftriangles, 0);
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            double err = (1. / 3.) * B[i].norm();
            Elem myelem = mesh.data.trianglelist[i];
            // make sure the vertices are in counterclockwise order
            if (Vector::doubleOrientedArea(mesh.data.pointlist[myelem[0]],
                                           mesh.data.pointlist[myelem[1]],
                                           mesh.data.pointlist[myelem[2]]) < 0) {
                std::swap(myelem[0], myelem[1]);
            }
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t id1 = myelem[j];
                uint32_t id2 = myelem[(j + 1) % 3];
                Vector v1 = mesh.data.pointlist[id1];
                Vector v2 = mesh.data.pointlist[id2];
                double len = Vector::distance(v1, v2);
                Vector n = (v2 - v1).normal().normalize();
                uint64_t edge = id1 < id2 ? (uint64_t)id1 << 32 | id2 : (uint64_t)id2 << 32 | id1;
                uint32_t adjid1 = edge_to_triangles[edge].elems[0];
                uint32_t adjid2 = edge_to_triangles[edge].elems[1];
                Vector Bp = B[adjid1 == i ? adjid2 : adjid1];
                err += (n ^ (Bp - B[i])) * len;
            }
            elem_field_errors[i] = err;
        }

        printf("computing median field error\n");
        // compute min and max field error
        double min_err = std::numeric_limits<double>::max();
        double max_err = std::numeric_limits<double>::min();
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            min_err = std::min(min_err, elem_field_errors[i]);
            max_err = std::max(max_err, elem_field_errors[i]);
        }

        // square root of 2 / 2
        Vector myS = Vector(1, 0).normalize();

        printf("building fem matrix\n");
        // surface integral
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            for (uint8_t j = 0; j < adjelems_count[i]; j++) {
                uint32_t v1, v2, v3;
                v1 = i;
                Elem myelem = mesh.data.trianglelist[adjelems_ids[i][j]];
                Vector myB = B[adjelems_ids[i][j]];
                if (i == mesh.data.trianglelist[adjelems_ids[i][j]][0]) {
                    v2 = myelem[1];
                    v3 = myelem[2];
                } else if (i == myelem[1]) {
                    v2 = myelem[2];
                    v3 = myelem[0];
                } else if (i == myelem[2]) {
                    v2 = myelem[0];
                    v3 = myelem[1];
                } else {
                    nexit("error: vertex not found in element");
                }

                double oriented_area = mesh.data.getDoubleOrientedArea(v1, v2, v3);
                
                if (oriented_area < 0) {
                    std::swap(v2, v3);
                }

                double area = mesh.data.getDoubleOrientedArea(v1, v2, v3);

                // elements are only added if they are in the upper triangle because the matrix is symmetric and this saves half the memory
                double b1 = (mesh.data.pointlist[v2].y - mesh.data.pointlist[v3].y) / area;
                double b2 = (mesh.data.pointlist[v3].y - mesh.data.pointlist[v1].y) / area;
                double b3 = (mesh.data.pointlist[v1].y - mesh.data.pointlist[v2].y) / area;
                double c1 = (mesh.data.pointlist[v3].x - mesh.data.pointlist[v2].x) / area;
                double c2 = (mesh.data.pointlist[v1].x - mesh.data.pointlist[v3].x) / area;
                double c3 = (mesh.data.pointlist[v2].x - mesh.data.pointlist[v1].x) / area;
                // coo.add_elem(i, v1, (area * (b1 * b1 + c1 * c1)) / (2 * adjelems_props[i][j].mu));
                // coo.add_elem(i, v2, (area * (b2 * b1 + c2 * c1)) / (2 * adjelems_props[i][j].mu));
                // coo.add_elem(i, v3, (area * (b3 * b1 + c3 * c1)) / (2 * adjelems_props[i][j].mu));

                double K1 = - 0.5 * myS.x * myB.x * myB.x + 0.5 * myS.x * myB.y * myB.y - myS.y * myB.x * myB.y;
                double K2 = - 0.5 * myS.y * myB.y * myB.y + 0.5 * myS.y * myB.x * myB.x - myS.x * myB.x * myB.y;
                // double K1 = 1;
                // double K2 = 1;

                double err = elem_field_errors[adjelems_ids[i][j]];
                
                double Wi = limit(map(sqrt(fabs(err)), min_err, max_err, 1, 1000), 1, 1000);
                elem_field_weights[adjelems_ids[i][j]] = Wi;

                // if (v1 >= i) coo(i, v1) += (Wi * Wi * area * (K1 * b1 + K2 * c1) * (K1 * b1 + K2 * c1));
                // if (v2 >= i) coo(i, v2) += (Wi * Wi * area * (K1 * b2 + K2 * c2) * (K1 * b1 + K2 * c1));
                // if (v3 >= i) coo(i, v3) += (Wi * Wi * area * (K1 * b3 + K2 * c3) * (K1 * b1 + K2 * c1));

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 * K1 + c1 * c1 * K2 + 4 * b1 * c1 * K1 * K2);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b2 * b1 * K1 + c2 * c1 * K2 + 2 * K1 * K2 * (b1 * c2 + b2 * c1));
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b3 * b1 * K1 + c3 * c1 * K2 + 2 * K1 * K2 * (b1 * c3 + b3 * c1));

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 * K1 + c1 * c1 * K2 + 2 * K1 * K2 * b1 * c1);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b2 * b1 * K1 + c2 * c1 * K2 + 2 * K1 * K2 * b2 * c1);
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b3 * b1 * K1 + c3 * c1 * K2 + 2 * K1 * K2 * b3 * c1);

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 * K1 + c1 * c1 * K2 + 2 * K1 * K2 * c1 * b1);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b2 * b1 * K1 + c2 * c1 * K2 + 2 * K1 * K2 * c2 * b1);
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b3 * b1 * K1 + c3 * c1 * K2 + 2 * K1 * K2 * c3 * b1);

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 * K1 + c1 * c1 * K2);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b2 * b1 * K1 + c2 * c1 * K2);
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b3 * b1 * K1 + c3 * c1 * K2);
                
                if (v1 >= i) coo(i, v1) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b1 + K2 * c1);
                if (v2 >= i) coo(i, v2) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b2 + K2 * c2);
                if (v3 >= i) coo(i, v3) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b3 + K2 * c3);

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 + c1 * c1);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b1 * b2 + c1 * c2);
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b1 * b3 + c1 * c3);
            }
        }

        printf("dirichlet boundary conditions\n");
        // std::vector<double> test(mesh.data.numberofpoints, 0);
        // // we have to set a 1 dirichlet boundary condition for all the vertices inside the integration region
        // for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
        //     Vector mypoint = mesh.data.pointlist[i];
        //     bool inside_integration_region = integration_region.contains(mypoint, true, mesh.epsilon);
        //     bool inside_boundary_region = boundary_region.contains(mypoint, true, mesh.epsilon);
        //     bool inside_boundary_region_hole = false;
        //     for (auto& hole : boundary_region_holes) {
        //         if (hole.contains(mypoint, false, mesh.epsilon)) {
        //             inside_boundary_region_hole = true;
        //             break;
        //         }
        //     }
        //     // OPTIMIZE ---------------------------------------------------------------------------------------------
        //     if (inside_integration_region) {
        //         for (auto& e : coo.elems) {
        //             uint32_t col = e.first >> 32;
        //             uint32_t row = e.first & 0xFFFFFFFF;
        //             double old_value = e.second;
        //             // diagonal
        //             if (row == i && col == i) {
        //                 e.second = 1;
        //                 b[i] = 1;
        //                 //b_dirichlet_mask[i] = true;
        //                 continue;
        //             }
        //             // row < col the matrix is stored in the upper triangular part
        //             if (col == i) {
        //                 if (!b_dirichlet_mask[i]) b[row] -= old_value;
        //                 e.second = 0;
        //             }
        //             if (row == i) {
        //                 e.second = 0;
        //             }
        //             std::swap(row, col);
        //             if (col == i) {
        //                 if (!b_dirichlet_mask[i]) b[col] -= old_value;
        //             }
        //         }
        //         test [i] = 1;
        //     } else if (!inside_boundary_region || inside_boundary_region_hole) {
        //         for (auto& e : coo.elems) {
        //             uint32_t row = e.first >> 32;
        //             uint32_t col = e.first & 0xFFFFFFFF;
        //             if (row == i || col == i) {
        //                 e.second = 0;
        //             }
        //             if (row == i && col == i) {
        //                 e.second = 1;
        //             }
        //         }
        //         test [i] = 2;
        //     }
        //     if (i % 1000 == 0) {
        //         printf("i: %d\n", i);
        //     }
        // }

        for (auto& [id, value] : coo.elems) {
            
        }

        MatCSRSymmetric Ai(coo);
        auto g = CV(mesh.data.numberofpoints);

        preconditionedSSORConjugateGradientSolver(Ai, b, g, 1.5, 1e-8, 10000);

        NodeScalarPlotToFile(10000, 10000, g.val, "g.png");
        // ElemScalarPlotToFile(10000, 10000, elem_field_errors, "errors.png");
        // NodeScalarPlotToFile(10000, 10000, test, "test.png");        // A.plot("A.png");
        // ElemScalarPlotToFile(10000, 10000, elem_field_weights, "weights.png");

        // compute grad(g)
        std::vector<Vector> nabla_g(mesh.data.numberoftriangles, {0, 0});
        mesh.computeGrad(nabla_g, g);

        std::vector<double> nabla_g_norm(mesh.data.numberoftriangles, 0);
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            nabla_g_norm[i] = nabla_g[i].norm();
        }

        ElemScalarPlotToFile(10000, 10000, nabla_g_norm, "nabla_g.png");
        // ElemScalarPlotToFile(10000, 10000, boundary_elements_mask, "boundary.png");

        // integrate force
        Vector force = {0, 0};

        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            Vector vertices[3] = {
                mesh.data.pointlist[mesh.data.trianglelist[i][0]],
                mesh.data.pointlist[mesh.data.trianglelist[i][1]],
                mesh.data.pointlist[mesh.data.trianglelist[i][2]]
            };

            bool inside = true;
            for (uint32_t j = 0; j < 3; j++) {
                bool in_integration_region = integration_region.contains(vertices[j], false, mesh.epsilon);
                bool in_boundary_region = boundary_region.contains(vertices[j], true, mesh.epsilon);
                bool in_boundary_region_hole = false;
                for (auto& hole : boundary_region_holes) {
                    if (hole.contains(vertices[j], false, mesh.epsilon)) {
                        in_boundary_region_hole = true;
                        break;
                    }
                }

                if (in_integration_region || !in_boundary_region || in_boundary_region_hole) {
                    inside = false;
                    break;
                }
            }

            if (inside) {
                uint32_t v1 = mesh.data.trianglelist[i][0];
                uint32_t v2 = mesh.data.trianglelist[i][1];
                uint32_t v3 = mesh.data.trianglelist[i][2];

                Vector p1 = mesh.data.pointlist[v1];
                Vector p2 = mesh.data.pointlist[v2];
                Vector p3 = mesh.data.pointlist[v3];

                Vector my_nabla_g = nabla_g[i];
                Vector my_B = B[i];

                double area = Vector::area(p1, p2, p3);

                force += (my_nabla_g * 0.5 * my_B.normSquared() - my_B * (my_nabla_g * my_B)) * area;
            }
        }
        force *= 1 / MU_0;
        return force;

        // EGGSHELL
        /*
        // debug plot
        // create the image
        uint32_t width = 10000;
        uint32_t height = 10000;
        bool plotMesh = true;
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // get mesh enclosing rectangle
        float min_x = -1.1 * mesh.radius;
        float min_y = -1.1 * mesh.radius;
        float max_x = 1.1 * mesh.radius;
        float max_y = 1.1 * mesh.radius;

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

        // render
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Vector v1 = mesh.data.pointlist[e[0]];
            Vector v2 = mesh.data.pointlist[e[1]];
            Vector v3 = mesh.data.pointlist[e[2]];

            if (Vector::distance(v1, Vector(0, 0)) > mesh.radius + mesh.epsilon ||
                Vector::distance(v2, Vector(0, 0)) > mesh.radius + mesh.epsilon || 
                Vector::distance(v3, Vector(0, 0)) > mesh.radius + mesh.epsilon) {
                continue;
            }

            // draw mesh edges
            // draw lines from v1 to v2
            if (plotMesh) {
                cv::line(image, cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                                cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                                cv::Scalar(255, 255, 255), 1);
                cv::line(image, cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                                cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset),
                                cv::Scalar(255, 255, 255), 1);
                cv::line(image, cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset),
                                cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                                cv::Scalar(255, 255, 255), 1);
            }

            // level curves
            
        }

        // draw the geometry
        // draw the segments
        for (DrawingSegment s : mesh.drawing.segments) {
            cv::line(image, cv::Point(x_scale * mesh.drawing.points[s.p1].x + x_offset,
                                       y_scale * mesh.drawing.points[s.p1].y + y_offset), 
                             cv::Point(x_scale * mesh.drawing.points[s.p2].x + x_offset, 
                                       y_scale * mesh.drawing.points[s.p2].y + y_offset),
                             cv::Scalar(255, 255, 255), 1);
        }

        // plot red x at the center of the integration region
        Vector icenter = {0, 0};
        for (auto point : integration_region.points) {
            icenter += point;
        }
        icenter /= integration_region.points.size();
        cv::line(image, cv::Point(x_scale * icenter.x + x_offset - 5, y_scale * icenter.y + y_offset - 5),
                        cv::Point(x_scale * icenter.x + x_offset + 5, y_scale * icenter.y + y_offset + 5),
                        cv::Scalar(0, 0, 255), 1);
        cv::line(image, cv::Point(x_scale * icenter.x + x_offset - 5, y_scale * icenter.y + y_offset + 5),
                        cv::Point(x_scale * icenter.x + x_offset + 5, y_scale * icenter.y + y_offset - 5),
                        cv::Scalar(0, 0, 255), 1);

        // find all the elements that are on the exterior of the boundary of the polygon
        // get list of polygon segments and their normals
        struct SegmentNormal {
            uint32_t p1;
            uint32_t p2;
            Vector normal;
        };
        std::vector<SegmentNormal> segment_normals;
        for (uint32_t i = 0; i < integration_region.points.size(); i++) {
            uint32_t id1 = i;
            uint32_t id2 = (i + 1) % integration_region.points.size();
            Vector p1 = integration_region.points[id1];
            Vector p2 = integration_region.points[id2];
            Vector n = (p2 - p1).normal().normalize();
            // we need the outward normal
            Vector test = Vector::midPoint(p1, p2) + (n * mesh.epsilon * 0.1);  // 0.1 just to make sure that we get no false negatives
            if (integration_region.contains(test)) {
                n = n * -1;
            }
            segment_normals.push_back(
                {id1, id2, n}
            );

            // plot the normals as arrows
            Vector mid = Vector::midPoint(p1, p2);
            cv::arrowedLine(image, cv::Point(x_scale * mid.x + x_offset, y_scale * mid.y + y_offset),
                                   cv::Point(x_scale * (mid.x + n.x * 0.1) + x_offset, y_scale * (mid.y + n.y * 0.1) + y_offset),
                                   cv::Scalar(0, 255, 0), 1);
        }
        
        // find all the nodes that are on the boundary of the polygon
        std::unordered_map<uint32_t, Vector> bnodes;
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector myp = mesh.data.pointlist[i];
            for (auto segment_normal : segment_normals) {
                Vector p1 = integration_region.points[segment_normal.p1];
                Vector p2 = integration_region.points[segment_normal.p2];
                double dist = Segment::pointSegmentDistance(myp, p1, p2);
                if (dist < mesh.epsilon * 0.1) {
                    bnodes[i] = (bnodes[i] + segment_normal.normal) / 2;
                }
            }
        }

        // plot bnodes as blue arrows
        for (auto bnode : bnodes) {
            uint32_t id = bnode.first;
            Vector myp = mesh.data.pointlist[id];
            Vector n = bnode.second;
            cv::arrowedLine(image, cv::Point(x_scale * myp.x + x_offset, y_scale * myp.y + y_offset),
                                   cv::Point(x_scale * (myp.x + n.x * 0.1) + x_offset, y_scale * (myp.y + n.y * 0.1) + y_offset),
                                   cv::Scalar(255, 0, 0), 1);
        }

        // find all the elements that are on the boundary of the polygon
        struct BElem {
            uint32_t id;
            Vector normal;
            double area;
            double length;
            Vector B;
        };

        std::vector<BElem> belems;
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            Elem e = mesh.data.trianglelist[i];
            uint32_t common_nodes[2] = {0, 0};
            uint32_t common_nodes_count = 0;
            for (auto bnode : bnodes) {
                uint32_t id = bnode.first;
                Vector myp = mesh.data.pointlist[id];
                Vector n = bnode.second;
                if (e[0] == id || e[1] == id || e[2] == id) {
                    common_nodes[common_nodes_count] = id;
                    common_nodes_count++;
                    if (common_nodes_count > 3) {
                        nexit("too many common nodes, this should never happen");
                    }
                }
            }
            switch (common_nodes_count) {
                case 0:
                    break;
                case 3:  // this triangle has one node in the intersection of two segments of the polygon and is inside the polygon
                    break;
                case 1: {
                    // check both other nodes are outside the polygon
                    Vector onode = mesh.data.pointlist[common_nodes[0]];
                    Vector anode, bnode;
                    if (e[0] == common_nodes[0]) {
                        anode = mesh.data.pointlist[e[1]];
                        bnode = mesh.data.pointlist[e[2]];
                    } else if (e[1] == common_nodes[0]) {
                        anode = mesh.data.pointlist[e[0]];
                        bnode = mesh.data.pointlist[e[2]];
                    } else if (e[2] == common_nodes[0]) {
                        anode = mesh.data.pointlist[e[0]];
                        bnode = mesh.data.pointlist[e[1]];
                    }
                    Vector n = bnodes[common_nodes[0]];
                    if (!integration_region.contains(anode) && !integration_region.contains(bnode)) {
                        // we have a boundary element
                        double area = fabs(Vector::doubleOrientedArea(anode, onode, bnode) * 0.5);
                        double length = Vector::distance(anode, bnode);
                        belems.push_back(
                            {i, n, area, -length, B[i]} // negative length because of parametric integration
                        );
                        // plot the boundary element as a green filled triangle
                        cv::Point points[1][3];
                        points[0][0] = cv::Point(x_scale * anode.x + x_offset, y_scale * anode.y + y_offset);
                        points[0][1] = cv::Point(x_scale * onode.x + x_offset, y_scale * onode.y + y_offset);
                        points[0][2] = cv::Point(x_scale * bnode.x + x_offset, y_scale * bnode.y + y_offset);
                        const cv::Point* ppt[1] = {points[0]};
                        int npt[] = {3};
                        cv::fillPoly(image, ppt, npt, 1, cv::Scalar(0, 255, 0));
                    }

                    break;
                }
                case 2: {
                    // check the other node is outside the polygon
                    Vector anode = mesh.data.pointlist[common_nodes[0]];
                    Vector bnode = mesh.data.pointlist[common_nodes[1]];
                    Vector onode;
                    if (e[0] != common_nodes[0] && e[0] != common_nodes[1]) {
                        onode = mesh.data.pointlist[e[0]];
                    } else if (e[1] != common_nodes[0] && e[1] != common_nodes[1]) {
                        onode = mesh.data.pointlist[e[1]];
                    } else if (e[2] != common_nodes[0] && e[2] != common_nodes[1]) {
                        onode = mesh.data.pointlist[e[2]];
                    }
                    if (!integration_region.contains(onode)) {
                        // we have a boundary element
                        double area = fabs(Vector::doubleOrientedArea(anode, onode, bnode) * 0.5);
                        double length = Vector::distance(anode, bnode);
                        belems.push_back(
                            {i, (bnodes[common_nodes[0]] + bnodes[common_nodes[1]]) / 2, area, length, B[i]} // positive length because of parametric integration
                        );
                        // plot the boundary element as a red filled triangle
                        cv::Point points[1][3];
                        points[0][0] = cv::Point(x_scale * anode.x + x_offset, y_scale * anode.y + y_offset);
                        points[0][1] = cv::Point(x_scale * onode.x + x_offset, y_scale * onode.y + y_offset);
                        points[0][2] = cv::Point(x_scale * bnode.x + x_offset, y_scale * bnode.y + y_offset);
                        const cv::Point* ppt[1] = {points[0]};
                        int npt[] = {3};
                        cv::fillPoly(image, ppt, npt, 1, cv::Scalar(0, 0, 255));
                    }

                    break;
                }
                default:
                    nexit("too many common nodes, this should never happen");
            }
        }

        // compute the force integral
        Vector force_integral = Vector(0, 0);

        for (auto belem : belems) {
            Elem e = mesh.data.trianglelist[belem.id];
            Vector n = belem.normal;
            double area = belem.area;
            double length = belem.length;
            Vector p1 = mesh.data.pointlist[e[0]];
            Vector p2 = mesh.data.pointlist[e[1]];
            Vector p3 = mesh.data.pointlist[e[2]];
            Vector Bm = belem.B;
            Vector sigma = {
                (+ 0.5 * (Bm.x * Bm.x - Bm.y * Bm.y) * n.x + Bm.x * Bm.y * n.y) * (1 / MU_0),
                (- 0.5 * (Bm.x * Bm.x - Bm.y * Bm.y) * n.y + Bm.x * Bm.y * n.x) * (1 / MU_0)
            };
            force_integral += sigma * area;
            // force_integral += sigma * length * 0.5;
        }

        // save the image
        cv::imwrite("kek.png", image);

        return force_integral;
        */
    }
}