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
#include "../algebra/coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/solvers.hpp"
#include "../algebra/math.hpp"

namespace nikfemm {
#ifdef NIKFEMM_USE_OPENCV
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

        std::vector<float> A_sorted(A.size());
        for (uint32_t i = 0; i < A.size(); i++) {
            A_sorted[i] = A[i];
        }
        std::sort(A_sorted.begin(), A_sorted.end());
        float max_A = A_sorted[0.9 * A_sorted.size()];
        float min_A = A_sorted[0.1 * A_sorted.size()];

        nloginfo("max_A: %f, min_A: %f", max_A, min_A);
        

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

        // flip image horizontally
        cv::flip(image, image, 0);

        cv::imshow("A", image);
        nloginfo("showing image");
        cv::waitKey(0);
        nloginfo("done");
    }

    void MagnetostaticSimulation::AplotToFile(uint32_t width, uint32_t height, std::string filename) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        if (!image.data) {
            nexit("Could not create image");
        }

        AplotRend(&image, width, height);

        // flip image horizontally
        cv::flip(image, image, 0);

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

        nloginfo("max_A: %f, min_A: %f", max_A, min_A);
        

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

        // flip image horizontally
        cv::flip(image, image, 0);

        cv::imshow("A", image);
        nloginfo("showing image");
        cv::waitKey(0);
        nloginfo("done");
        }

    void MagnetostaticSimulation::NodeScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh, bool plotRegions) {
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        if (!image.data) {
            nexit("Could not create image");
        }

        NodeScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        // flip image horizontally
        cv::flip(image, image, 0);

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

        nloginfo("max B: %f", max_B);
        nloginfo("min B: %f", min_B);
        

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

    void MagnetostaticSimulation::Bplot(uint32_t width, uint32_t height, bool plotMesh, bool plotRegions, double maxB, double minB, bool waitkey, bool plotCurves, uint32_t curve_number) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        BplotRend(&image, width, height, plotMesh, plotRegions, maxB, minB, plotCurves, curve_number);

        // flip image horizontally
        cv::flip(image, image, 0);

        // show the image
        cv::imshow("B", image);
        // continue if image is closed
        if (waitkey) {
            cv::waitKey(0);
        } else {
            cv::waitKey(1);
        }
    }

    void MagnetostaticSimulation::BplotToFile(uint32_t width, uint32_t height, std::string filename, bool plotMesh, bool plotRegions, double maxB, double minB, bool plotCurves, uint32_t curve_number) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        BplotRend(&image, width, height, plotMesh, plotRegions, maxB, minB, plotCurves, curve_number);

        // flip image horizontally
        cv::flip(image, image, 0);

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

        nloginfo("max_Scalar: %f, min_Scalar: %f", max_Scalar, min_Scalar);
        

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

        // flip image horizontally
        cv::flip(image, image, 0);

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

        // flip image horizontally
        cv::flip(image, image, 0);

        // save the image
        cv::imwrite(filename, image);
    }
#endif

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
                mu[i] = (props[i]->getMu(Bmag) * kalman_scemo) + ((mu[i] + magnetostatic_materials::vacuum * residual * residual) * (1 - kalman_scemo)); 
                // mu[i] += (props[i]->getMu(Bmag) - mu[i]) * 0.1;  // mu += (mu_new - mu) * 0.1
                // mu[i] += props[i]->getMu(Bmag);  // pure newton-raphson
            }
        }
    }

    System<MagnetostaticNonLinearExpression> MagnetostaticSimulation::generateSystem(bool refine, double max_triangle_area, int min_angle) {
        // get time in milliseconds

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        // Circle smallest_circle = Circle::getMinimumEnclosingCircle(mesh.drawing.points);
        Circle smallest_circle = Circle::getEnclosingCircle(mesh.drawing.points);
        if (mesh.drawing.points.size() == 0) {
            smallest_circle.radius = 1;
        }
        double coeff = 1.5;  // make the circle [coeff] times the size of the smallest circle
        // set simulation offset and boundary radius
        mesh.center = smallest_circle.center;
        mesh.radius = smallest_circle.radius * coeff;
        // translate everything to the origin
        mesh.drawing.translate(-mesh.center);
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(Vector(0, 0), smallest_circle.radius * coeff);
        double circumferential_length = boundary_circle.circumference();
        uint32_t boundary_points = (uint32_t)((circumferential_length / sqrt((max_triangle_area * 4) / sqrt(3))) * 2);
        nloginfo("boundary points: %d", boundary_points);
        mesh.drawing.drawCircle(boundary_circle, boundary_points);
        // add region near the edge of the circle
        mesh.drawing.drawRegion(Vector(boundary_circle.radius * 0.9, 0), {0, {0, 0}, magnetostatic_materials::air});
        // add the boundary 
        // mesh.drawing.plot();
        #ifdef NIK_REFINE_MAGNETS
            mesh.refineMeshAroundMagnets();
        #endif

        if (refine) {
            mesh.drawing.addRefiningPoints();
        }

        mesh.mesh(max_triangle_area, min_angle);
        // mesh.plot();
        mesh.addKelvinBoundaryConditions(boundary_points);
        mesh.computeEpsilon();
        // mesh.plot();
        nloginfo("the mesh has %u nodes and %u elements", mesh.data.numberofpoints, mesh.data.numberoftriangles);
        

        auto system = mesh.getFemSystem();
        mesh.addDirichletInfiniteBoundaryConditions(system);
        return system;
    }

    void MagnetostaticSimulation::solve(System<MagnetostaticNonLinearExpression>& system) {
        A = std::vector<double>(mesh.data.numberofpoints);
        B = std::vector<Vector>(mesh.data.numberoftriangles, {0, 0});
        std::vector<float> mu(mesh.data.numberoftriangles, 0);
        std::vector<const MagnetostaticProp*> props(mesh.data.numberoftriangles);

        // fill props
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            props[i] = mesh.drawing.getRegionPtrFromId(mesh.data.triangleattributelist[i]);
        }

        MagnetostaticMatCSRSymmetric FemMat(system.A);
        std::vector<double> b(system.b.size());

        // initialize mu
        for (uint32_t i = 0; i < B.size(); i++) {
            mu[i] = props[i]->getMu(0);
        }
        FemMat.updateFromMu(mu);
        for (uint32_t i = 0; i < system.b.size(); i++) {
            b[i] = system.b[i].evaluate(mu);
        }

        // check if magnetostatic_materials are all linear
        bool all_linear = true;
        for (uint32_t i = 0; i < props.size(); i++) {
            if (!props[i]->isLinear()) {
                all_linear = false;
                break;
            }
        }

        auto tstart = std::chrono::high_resolution_clock::now();
        if (all_linear) {
            nloginfo("all magnetostatic_materials are linear");
            preconditionedSSORConjugateGradientSolver(FemMat, b, A, 1.5, 1e-6, 100000);
            // preconditionedIncompleteCholeskyConjugateGradientSolver(FemMat, b, A, 1e-6, 100000);
        } else {
            nloginfo("nonlinear magnetostatic_materials detected, starting non linear newton solver");
            std::vector<double> r(system.b.size());  // residual

            double residual = 1e10;
            for (uint32_t i = 0; i < 500; i++) {
                // conjugateGradientSolver(FemMat, b, A, 1e-7, 10000);
                // preconditionedJacobiConjugateGradientSolver(FemMat, b, A, 1e-7, 1000);
                preconditionedSSORConjugateGradientSolver(FemMat, b, A, 1.5, 1e-7, 50);
                // preconditionedIncompleteCholeskyConjugateGradientSolver(FemMat, b, A, 1e-7, 1000);
                // mesh.Aplot(A);
                // mesh.Bplot(B);

                mesh.computeCurl(B, A);
                // mesh.Bplot(B);

                // check if the solution is correct
                MagnetostaticSimulation::updateMu(props, mu, B, residual, i);
                FemMat.updateFromMu(mu);
                for (uint32_t i = 0; i < system.b.size(); i++) {
                    b[i] = system.b[i].evaluate(mu);
                }
                mult(r, FemMat, A);
                sub(r, b, r);
                residual = norm(r);
                if (residual < 1e-7) {
                    nloginfo("Converged in %d iterations", i);
                    break;
                }
                // nloginfo("%.17g, %.17g; ", K, residual);
                nloginfo("nonlinear iteration %d, residual: %.17g", i, residual);
                // print mu
                // nloginfo("%f,", residual);
                fflush(stdout);
            }
        }
        auto tend = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count() / 1000.0);

        mesh.computeCurl(B, A);
    }

    Polygon MagnetostaticSimulation::getInnermostPolygon(Vector p) {
        // find the polygon that contains p
        Polygon integration_region;
        Vector translated_p = p - mesh.center;
        std::vector<Polygon> polygons_that_contain_p;
        for (auto polygon : mesh.drawing.polygons) {
            if (polygon.contains(translated_p)) {
                polygons_that_contain_p.push_back(polygon);
            }
        }

        // if there is more than one polygon that contains p then we have to find the polygon that isn't contained by any other polygon that contains p
        if (polygons_that_contain_p.size() == 0) {
            nexit("ERROR: no polygon contains p");
        } else if (polygons_that_contain_p.size() == 1) {
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

        return integration_region;
    }

    auto MagnetostaticSimulation::getSurroundingRegionBlockIntegralAssets(Vector p) {
        // find the polygon that contains p
        Polygon integration_region = getInnermostPolygon(p);
        Polygon boundary_region;
        std::vector<Polygon> boundary_region_holes;

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

        // COMPUTE FEM WEIGHTS
        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(mesh.data.numberofpoints);
        auto adjelems_count = std::vector<uint8_t>(mesh.data.numberofpoints, 0);

        nloginfo("computing adjelems_ids");
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = mesh.data.trianglelist[i][j];
                adjelems_ids[myid][adjelems_count[myid]++] = i;
            }
        }

        nloginfo("computing elemadjelems_ids");
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

        nloginfo("computing field error");
        // compute field error for each element
        std::vector<double> elem_field_errors = std::vector<double>(mesh.data.numberoftriangles, 0);
        // std::vector<double> elem_field_weights = std::vector<double>(mesh.data.numberoftriangles, 0);
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

        nloginfo("computing median field error");
        // compute min and max field error
        double min_err = std::numeric_limits<double>::max();
        double max_err = std::numeric_limits<double>::min();
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            min_err = std::min(min_err, elem_field_errors[i]);
            max_err = std::max(max_err, elem_field_errors[i]);
        }

        nloginfo("find position of vertices");
        std::vector<bool> vertex_inside_integration_region(mesh.data.numberofpoints);
        std::vector<bool> vertex_inside_integration_region_with_boundary(mesh.data.numberofpoints);
        std::vector<bool> vertex_inside_boundary_region(mesh.data.numberofpoints);
        std::vector<bool> vertex_inside_boundary_region_hole(mesh.data.numberofpoints);

        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector mypoint = mesh.data.pointlist[i];
            bool inside_integration_region_with_boundary = integration_region.contains(mypoint, true, mesh.epsilon);
            bool inside_integration_region = integration_region.contains(mypoint, false, mesh.epsilon);
            bool inside_boundary_region = boundary_region.contains(mypoint, true, mesh.epsilon);
            bool inside_boundary_region_hole = false;
            for (auto& hole : boundary_region_holes) {
                if (hole.contains(mypoint, false, mesh.epsilon)) {
                    inside_boundary_region_hole = true;
                    break;
                }
            }

            vertex_inside_integration_region[i] = inside_integration_region;
            vertex_inside_integration_region_with_boundary[i] = inside_integration_region_with_boundary;
            vertex_inside_boundary_region[i] = inside_boundary_region;
            vertex_inside_boundary_region_hole[i] = inside_boundary_region_hole;
        }

        return SurroundingRegionBlockIntegralAssets{
            adjelems_ids,
            adjelems_count, 
            elem_field_errors, 
            min_err, 
            max_err, 
            vertex_inside_integration_region,
            vertex_inside_integration_region_with_boundary,
            vertex_inside_boundary_region,
            vertex_inside_boundary_region_hole
        };
    }

    Vector MagnetostaticSimulation::computeForceIntegrals(SurroundingRegionBlockIntegralAssets assets, Vector p) {
        auto start = std::chrono::high_resolution_clock::now();
        
        auto [
            adjelems_ids, 
            adjelems_count, 
            elem_field_errors, 
            min_err, 
            max_err, 
            vertex_inside_integration_region,
            vertex_inside_integration_region_with_boundary,
            vertex_inside_boundary_region,
            vertex_inside_boundary_region_hole
        ] = assets;

        auto coo = MatCOOSymmetric<double>(mesh.data.numberofpoints);
        auto b = std::vector<double>(mesh.data.numberofpoints);
        auto b_dirichlet_mask = std::vector<bool>(mesh.data.numberofpoints, false);
        // since the stiffness matrix is symmetric, this function only computes the upper triangular part

        // square root of 2 / 2
        Vector myS = Vector(1, 0).normalize();

        nloginfo("building fem matrix");
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

                double K1 = - 0.5 * myS.x * myB.x * myB.x + 0.5 * myS.x * myB.y * myB.y - myS.y * myB.x * myB.y;
                double K2 = - 0.5 * myS.y * myB.y * myB.y + 0.5 * myS.y * myB.x * myB.x - myS.x * myB.x * myB.y;

                double err = elem_field_errors[adjelems_ids[i][j]];
                
                double Wi = limit(map(sqrt(fabs(err)), min_err, max_err, 1, 1000), 1, 1000);
                // elem_field_weights[adjelems_ids[i][j]] = Wi;
                
                if (v1 >= i) coo(i, v1) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b1 + K2 * c1);
                if (v2 >= i) coo(i, v2) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b2 + K2 * c2);
                if (v3 >= i) coo(i, v3) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b3 + K2 * c3);

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 + c1 * c1);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b1 * b2 + c1 * c2);
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b1 * b3 + c1 * c3);
            }
        }

        nloginfo("dirichlet boundary conditions");
        // std::vector<double> test(mesh.data.numberofpoints, 0);
        // we have to set a 1 dirichlet boundary condition for all the vertices inside the integration region

        uint32_t val = 0;
        for (auto& [id, value] : coo.elems) {
            uint32_t col = id >> 32;
            uint32_t row = id & 0xFFFFFFFF;

            bool col_in_zero_region = !vertex_inside_boundary_region[col] || vertex_inside_boundary_region_hole[col];
            bool col_in_one_region = vertex_inside_integration_region_with_boundary[col];

            double old_value = value;

            if (col_in_zero_region) {
                if (row == col) {
                    value = 1;
                    b[row] = 0;
                } else {
                    value = 0;
                }
            } else if (col_in_one_region) {
                if (row == col) {
                    value = 1;
                    b[row] = 1;
                } else {
                    value = 0;
                    b[row] -= old_value;
                }
            }
            val++;
        }

        MatCSRSymmetric Ai(coo);
        // Ai.write_to_file("Ad.txt");
        auto g = std::vector<double>(mesh.data.numberofpoints);

        // the error for this function varies a lot, so I just set the gradient to do at least 100 iterations and eventually stop at 1e-100
        preconditionedSSORConjugateGradientSolver(Ai, b, g, 1.5, 1e-100, 100);

        // NodeScalarPlotToFile(10000, 10000, g.val, "g.png");
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

        // ElemScalarPlotToFile(10000, 10000, nabla_g_norm, "nabla_g.png");
        // ElemScalarPlotToFile(10000, 10000, boundary_elements_mask, "boundary.png");

        // integrate force
        Vector force = {0, 0};

        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            bool inside = true;
            for (uint32_t j = 0; j < 3; j++) {
                uint32_t v = mesh.data.trianglelist[i][j];
                if (vertex_inside_integration_region[v] || !vertex_inside_boundary_region[v] || vertex_inside_boundary_region_hole[v]) {
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
        force *= 1 / (4 * PI * 1e-7);
        force *= depth;

        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("Time elapsed: %f ms", std::chrono::duration<double, std::milli>(end - start).count());
        return force;
    }

    double MagnetostaticSimulation::computeTorqueIntegral(SurroundingRegionBlockIntegralAssets assets, Vector p, Vector center) {
        auto start = std::chrono::high_resolution_clock::now();
        
        auto [
            adjelems_ids, 
            adjelems_count, 
            elem_field_errors, 
            min_err, 
            max_err, 
            vertex_inside_integration_region,
            vertex_inside_integration_region_with_boundary,
            vertex_inside_boundary_region,
            vertex_inside_boundary_region_hole
        ] = assets;

        auto coo = MatCOOSymmetric<double>(mesh.data.numberofpoints);
        auto b = std::vector<double>(mesh.data.numberofpoints);
        auto b_dirichlet_mask = std::vector<bool>(mesh.data.numberofpoints, false);
        // since the stiffness matrix is symmetric, this function only computes the upper triangular part

        nloginfo("building fem matrix");
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

                Vector barycenter = (mesh.data.pointlist[v1] + mesh.data.pointlist[v2] + mesh.data.pointlist[v3]) / 3;
                Vector r = barycenter - center;

                double K1 = - r.y * 0.5 * myB.normSquared() - myB.x * (r ^ myB);
                double K2 = r.x * 0.5 * myB.normSquared() - myB.y * (r ^ myB);

                double err = elem_field_errors[adjelems_ids[i][j]];
                
                double Wi = limit(map(sqrt(fabs(err)), min_err, max_err, 1, 1000), 1, 1000);
                // elem_field_weights[adjelems_ids[i][j]] = Wi;
                
                if (v1 >= i) coo(i, v1) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b1 + K2 * c1);
                if (v2 >= i) coo(i, v2) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b2 + K2 * c2);
                if (v3 >= i) coo(i, v3) += area * Wi * (K1 * b1 + K2 * c1) * (K1 * b3 + K2 * c3);

                // if (v1 >= i) coo(i, v1) += area * 0.5 * Wi * Wi * (b1 * b1 + c1 * c1);
                // if (v2 >= i) coo(i, v2) += area * 0.5 * Wi * Wi * (b1 * b2 + c1 * c2);
                // if (v3 >= i) coo(i, v3) += area * 0.5 * Wi * Wi * (b1 * b3 + c1 * c3);
            }
        }

        nloginfo("dirichlet boundary conditions");
        // std::vector<double> test(mesh.data.numberofpoints, 0);
        // we have to set a 1 dirichlet boundary condition for all the vertices inside the integration region

        uint32_t val = 0;
        for (auto& [id, value] : coo.elems) {
            uint32_t col = id >> 32;
            uint32_t row = id & 0xFFFFFFFF;

            bool col_in_zero_region = !vertex_inside_boundary_region[col] || vertex_inside_boundary_region_hole[col];
            bool col_in_one_region = vertex_inside_integration_region_with_boundary[col];

            double old_value = value;

            if (col_in_zero_region) {
                if (row == col) {
                    value = 1;
                    b[row] = 0;
                } else {
                    value = 0;
                }
            } else if (col_in_one_region) {
                if (row == col) {
                    value = 1;
                    b[row] = 1;
                } else {
                    value = 0;
                    b[row] -= old_value;
                }
            }
            val++;
        }

        MatCSRSymmetric Ai(coo);
        // Ai.write_to_file("Ad.txt");
        auto g = std::vector<double>(mesh.data.numberofpoints);

        // the error for this function varies a lot, so I just set the gradient to do at least 100 iterations and eventually stop at 1e-100
        preconditionedSSORConjugateGradientSolver(Ai, b, g, 1.5, 1e-100, 100);

        // NodeScalarPlotToFile(10000, 10000, g.val, "g_torque.png");
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

        // ElemScalarPlotToFile(10000, 10000, nabla_g_norm, "nabla_g_torque.png");
        // ElemScalarPlotToFile(10000, 10000, boundary_elements_mask, "boundary.png");

        // integrate force
        double torque = 0;

        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            bool inside = true;
            for (uint32_t j = 0; j < 3; j++) {
                uint32_t v = mesh.data.trianglelist[i][j];
                if (vertex_inside_integration_region[v] || !vertex_inside_boundary_region[v] || vertex_inside_boundary_region_hole[v]) {
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

                Vector barycenter = (p1 + p2 + p3) / 3;
                Vector r = barycenter - center;

                torque += ((r ^ my_nabla_g) * my_B.normSquared() * 0.5 - (r ^ my_B) * (my_B * my_nabla_g)) * area;
                // force += (my_nabla_g * 0.5 * my_B.normSquared() - my_B * (my_nabla_g * my_B)) * area;
            }
        }
        torque *= 1 / (4 * PI * 1e-7);
        torque *= depth;

        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("Time elapsed: %f ms", std::chrono::duration<double, std::milli>(end - start).count());
        return torque;
    }

    Vector MagnetostaticSimulation::computeForceIntegrals(Vector p) {
        auto assets = getSurroundingRegionBlockIntegralAssets(p);
        return computeForceIntegrals(assets, p);
    }

    double MagnetostaticSimulation::computeTorqueIntegral(Vector p, Vector center) {
        auto assets = getSurroundingRegionBlockIntegralAssets(p);
        return computeTorqueIntegral(assets, p, center);
    }

    StressTensor MagnetostaticSimulation::computeStressIntegral(Vector p, Vector center) {
        auto assets = getSurroundingRegionBlockIntegralAssets(p);
        Vector force = computeForceIntegrals(assets, p);
        double torque = computeTorqueIntegral(assets, p, center);
        return {force, torque};
    }

    double MagnetostaticSimulation::computeAreaIntegral(Vector p) {
        Polygon integration_region = getInnermostPolygon(p);

        double area = 0;
        
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            Vector p1 = mesh.data.pointlist[mesh.data.trianglelist[i][0]];
            Vector p2 = mesh.data.pointlist[mesh.data.trianglelist[i][1]];
            Vector p3 = mesh.data.pointlist[mesh.data.trianglelist[i][2]];

            if (integration_region.contains(p1) && integration_region.contains(p2) && integration_region.contains(p3)) {
                area += Vector::area(p1, p2, p3);
            }
        }

        return area;
    }

    double MagnetostaticSimulation::computeMass(Vector p, double density) {
        double area = computeAreaIntegral(p);
        return area * depth * density;
    }

    Vector MagnetostaticSimulation::computeBarycenter(Vector p) {
        Polygon integration_region = getInnermostPolygon(p);

        double area = 0;
        Vector barycenter = {0, 0};

        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            Vector p1 = mesh.data.pointlist[mesh.data.trianglelist[i][0]];
            Vector p2 = mesh.data.pointlist[mesh.data.trianglelist[i][1]];
            Vector p3 = mesh.data.pointlist[mesh.data.trianglelist[i][2]];

            if (integration_region.contains(p1) && integration_region.contains(p2) && integration_region.contains(p3)) {
                double triangle_area = Vector::area(p1, p2, p3);
                area += triangle_area;
                barycenter += ((p1 + p2 + p3) / 3) * triangle_area;
            }
        }

        return barycenter / area;
    }

    double MagnetostaticSimulation::computeInertiaMoment(Vector p, Vector center, double density) {
        Polygon integration_region = getInnermostPolygon(p);

        double moment = 0;

        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            Vector p1 = mesh.data.pointlist[mesh.data.trianglelist[i][0]];
            Vector p2 = mesh.data.pointlist[mesh.data.trianglelist[i][1]];
            Vector p3 = mesh.data.pointlist[mesh.data.trianglelist[i][2]];

            if (integration_region.contains(p1) && integration_region.contains(p2) && integration_region.contains(p3)) {
                double triangle_area = Vector::area(p1, p2, p3);
                Vector triangle_barycenter = (p1 + p2 + p3) / 3;
                Vector r = triangle_barycenter - center;
                moment += triangle_area * r.normSquared();
            }
        }

        return moment * depth * density;
    }
}