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

    void MagnetostaticSimulation::ScalarPlotRend(cv::Mat* image, double width, double height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
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
        // 99th percentile
        double max_Scalar = 0.00000000000000000000001;
        double min_Scalar = -0.00000000000000000000001;

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

    void MagnetostaticSimulation::ScalarPlot(uint32_t width, uint32_t height, std::vector<double>& scalar, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        ScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

        // show the image
        cv::imshow("B", image);
        // continue if image is closed
        cv::waitKey(0);
    }

    void MagnetostaticSimulation::ScalarPlotToFile(uint32_t width, uint32_t height, std::vector<double>& scalar, std::string filename, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        ScalarPlotRend(&image, width, height, scalar, plotMesh, plotRegions);

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

    MagnetostaticSystem MagnetostaticSimulation::generateSystem() {
        // get time in milliseconds
        mesh.drawing.addRefiningPoints();

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

        auto prop = mesh.drawing.getPolygonProp(integration_region);

        // compute simulation again
        MagnetostaticSimulation integral_simulation;
        integral_simulation.mesh.drawing.drawPolygon(integration_region);

        // add all points inside the polygon that contains p
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector my_point = mesh.data.pointlist[i];
            if (integration_region.contains(my_point)) {
                integral_simulation.mesh.drawing.drawPoint(my_point);
            }
        }

        // translate everything back
        integral_simulation.mesh.drawing.translate(mesh.center);
        integral_simulation.mesh.drawing.drawRegion(p, prop);
        // integral_simulation.mesh.drawing.plot(1000, 1000);

        // generate system of equations
        auto integral_system = integral_simulation.generateSystem();

        // for every point in the polygon that contains p set a dirichlet boundary condition so that the potential remains fixed
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector my_point = mesh.data.pointlist[i];
            if (integration_region.contains(my_point)) {
                // get potential from previous simulation
                double potential = 1;
                integral_simulation.mesh.addDirichletBoundaryConditions(integral_system, i, potential);
            }
        }

        MagnetostaticMatCSRSymmetric integral_FemMat(integral_system.coo);
        integral_FemMat.printCSR();

        // integral_system.b.write_to_file("bint.mtx");

        // solve system of equations
        integral_simulation.solve(integral_system);
        integral_simulation.BplotToFile(10000, 10000, "Bintegral.png", false, false);

        return Vector(0, 0);
    }
}