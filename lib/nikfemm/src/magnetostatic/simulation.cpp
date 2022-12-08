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
#include "../geometry/point.hpp"
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
            Point p = mesh.data.pointlist[i];

            if (Point::distance(p, Point(0, 0)) > mesh.radius + mesh.epsilon) {
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
                if (mesh.data.trianglelist[j].verts[0] == i || mesh.data.trianglelist[j].verts[1] == i || mesh.data.trianglelist[j].verts[2] == i) {
                    Point barycenter = {
                        (mesh.data.pointlist[mesh.data.trianglelist[j].verts[0]].x + mesh.data.pointlist[mesh.data.trianglelist[j].verts[1]].x + mesh.data.pointlist[mesh.data.trianglelist[j].verts[2]].x) / 3,
                        (mesh.data.pointlist[mesh.data.trianglelist[j].verts[0]].y + mesh.data.pointlist[mesh.data.trianglelist[j].verts[1]].y + mesh.data.pointlist[mesh.data.trianglelist[j].verts[2]].y) / 3
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

    void MagnetostaticSimulation::BplotRend(cv::Mat* image, double width, double height, bool plotMesh, bool plotRegions) {
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

        std::vector<double> B_mags;
        for (uint32_t i = 0; i < B.size(); i++) {
            B_mags.push_back(B[i].magnitude());
        }
        std::sort(B_mags.begin(), B_mags.end());
        double max_B = B_mags[0.9 * B_mags.size()];
        double min_B = B_mags[0.1 * B_mags.size()];

        printf("max B: %f\n", max_B);
        printf("min B: %f\n", min_B);

        // max_B = 1e-6;

        // render
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Point v1 = mesh.data.pointlist[e.verts[0]];
            Point v2 = mesh.data.pointlist[e.verts[1]];
            Point v3 = mesh.data.pointlist[e.verts[2]];

            if (Point::distance(v1, Point(0, 0)) > mesh.radius + mesh.epsilon ||
                Point::distance(v2, Point(0, 0)) > mesh.radius + mesh.epsilon || 
                Point::distance(v3, Point(0, 0)) > mesh.radius + mesh.epsilon) {
                continue;
            }

            // get the B vectors
            Vector Bv = B[i];
            cv::Scalar c = val2jet(Bv.magnitude(), min_B, max_B);
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
                             cv::Scalar(0, 0, 0), 1);
        }

        if (plotRegions) {
            // draw the regions
            for (DrawingRegion r : mesh.drawing.regions) {
                std::vector<cv::Point> points;
                
                Point pos = r.first;
                // draw white cross
                cv::line(*image, cv::Point(x_scale * pos.x + x_offset - 5, y_scale * pos.y + y_offset),
                                 cv::Point(x_scale * pos.x + x_offset + 5, y_scale * pos.y + y_offset),
                                 cv::Scalar(255, 255, 255), 1);
                cv::line(*image, cv::Point(x_scale * pos.x + x_offset, y_scale * pos.y + y_offset - 5),
                                 cv::Point(x_scale * pos.x + x_offset, y_scale * pos.y + y_offset + 5),
                                 cv::Scalar(255, 255, 255), 1);
                // draw an arrow for each magnet region
                /*
                if (mesh.drawing.region_map[r.second].M != Vector(0, 0)) {
                    Point start = r.first;
                    Point end = start + mesh.drawing.region_map[r.second].M;
                    cv::arrowedLine(*image, cv::Point(x_scale * start.x + x_offset, y_scale * start.y + y_offset),
                                            cv::Point(x_scale * end.x + x_offset, y_scale * end.y + y_offset),
                                            cv::Scalar(255, 255, 255), 1);
                }
                */
            }

            // draw all drawing points
            for (Point p : mesh.drawing.points) {
                cv::circle(*image, cv::Point(x_scale * p.x + x_offset, y_scale * p.y + y_offset), 2, cv::Scalar(255, 255, 255), -1);
            }
        }
    }

    void MagnetostaticSimulation::Bplot(uint32_t width, uint32_t height, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        BplotRend(&image, width, height, plotMesh, plotRegions);

        // show the image
        cv::imshow("B", image);
        // continue if image is closed
        cv::waitKey(0);
    }

    void MagnetostaticSimulation::BplotToFile(uint32_t width, uint32_t height, std::string filename, bool plotMesh, bool plotRegions) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        BplotRend(&image, width, height, plotMesh, plotRegions);

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
        double max_Scalar = 0.000000000000000000000000000000001;
        double min_Scalar = 0;

        printf("max Scalar: %.17g\n", max_Scalar);
        printf("min Scalar: %.17g\n", min_Scalar);

        // max_B = 1e-6;

        // render
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            // get the triangle
            Elem e = mesh.data.trianglelist[i];
            // get the vertices
            Point v1 = mesh.data.pointlist[e.verts[0]];
            Point v2 = mesh.data.pointlist[e.verts[1]];
            Point v3 = mesh.data.pointlist[e.verts[2]];

            if (Point::distance(v1, Point(0, 0)) > mesh.radius + mesh.epsilon ||
                Point::distance(v2, Point(0, 0)) > mesh.radius + mesh.epsilon || 
                Point::distance(v3, Point(0, 0)) > mesh.radius + mesh.epsilon) {
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
                
                Point pos = r.first;
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

    void MagnetostaticSimulation::solve() {
        // get time in milliseconds
        omp_set_num_threads(12);
        std::cout << "Number of threads in the current parallel region is " << omp_get_num_threads() << std::endl;
        mesh.drawing.addRefiningPoints();

        /* auto boundary */
        // find smallest enclosing circle using Welzl's algorithm
        auto start1 = std::chrono::high_resolution_clock::now();
        Circle smallest_circle = Circle::getMinimumEnclosingCircle(mesh.drawing.points);
        if (mesh.drawing.points.size() == 0) {
            smallest_circle.radius = 1;
        }
        // translate everything to the origin
        for (uint32_t i = 0; i < mesh.drawing.points.size(); i++) {
            mesh.drawing.points[i] = Point(mesh.drawing.points[i].x - smallest_circle.center.x, mesh.drawing.points[i].y - smallest_circle.center.y);
        }
        for (uint32_t i = 0; i < mesh.drawing.polygons.size(); i++) {
            for (uint32_t j = 0; j < mesh.drawing.polygons[i].points.size(); j++) {
                mesh.drawing.polygons[i].points[j] = Point(mesh.drawing.polygons[i].points[j].x - smallest_circle.center.x, mesh.drawing.polygons[i].points[j].y - smallest_circle.center.y);
            }
        }
        std::vector<DrawingRegion> translated_regions;
        for (DrawingRegion region : mesh.drawing.regions) {
            translated_regions.push_back(DrawingRegion(Point(region.first.x - smallest_circle.center.x, region.first.y - smallest_circle.center.y), region.second));
        }
        mesh.drawing.regions = translated_regions;
        // set simulation offset and boundary radius
        mesh.center = smallest_circle.center;
        mesh.radius = smallest_circle.radius * 2;
        // make circle double the size of the smallest circle
        Circle boundary_circle = Circle(Point(0, 0), 2 * smallest_circle.radius);
        double circumferential_length = boundary_circle.circumference();
        uint32_t boundary_points = (uint32_t)((circumferential_length / sqrt((MAX_TRIANGLE_AREA * 4) / sqrt(3))) * 2);
        printf("boundary points: %d\n", boundary_points);
        mesh.drawing.drawCircle(boundary_circle, boundary_points);
        // add region near the edge of the circle
        mesh.drawing.drawRegion(Point(boundary_circle.radius * 0.9, 0), {0, {0, 0}, materials::air});
        // add the boundary 
        // mesh.drawing.plot();
        mesh.refineMeshAroundMagnets();
        auto start2 = std::chrono::high_resolution_clock::now();
        mesh.mesh();
        auto start3 = std::chrono::high_resolution_clock::now();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        mesh.addKelvinBoundaryConditions(boundary_points);
        mesh.computeEpsilon();
        auto start4 = std::chrono::high_resolution_clock::now();
        #ifdef DEBUG_PRINT
        // mesh.plot();
        #endif
        printf("the mesh has %u nodes and %u elements\n", mesh.data.numberofpoints, mesh.data.numberoftriangles);
        BuildMatCOO<MagnetostaticNonLinearExpression> coo(mesh.data.numberofpoints);
        CV b(mesh.data.numberofpoints);
        A = CV(mesh.data.numberofpoints);
        B = std::vector<Vector>(mesh.data.numberoftriangles, {0, 0});
        std::vector<float> mu(mesh.data.numberoftriangles, 0);
        std::vector<const MagnetostaticProp*> props(mesh.data.numberoftriangles);

        // fill props
        for (uint32_t i = 0; i < mesh.data.numberoftriangles; i++) {
            props[i] = mesh.drawing.getRegionPtrFromId(mesh.data.triangleattributelist[i]);
        }

        auto start5 = std::chrono::high_resolution_clock::now();

        std::vector<double> Jm = mesh.getFemSystem(coo, b);
        ScalarPlotToFile(10000, 10000, Jm, "J.png", true, true);
        auto start6 = std::chrono::high_resolution_clock::now();
        mesh.addDirichletBoundaryConditions(coo, b);
        auto start7 = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        printf("coo matrix m: %lu, n: %lu, elems: %lu\n", coo.m, coo.n, coo.elems.size());
        #endif

        MagnetostaticMatCSRSymmetric FemMat(coo);
        // b.write_to_file("b");
        auto start8 = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        // csr.print();
        // coo.plot();
        // x.print();
        // b.print();
        #endif

        // initialize mu
        for (uint32_t i = 0; i < B.size(); i++) {
            mu[i] = props[i]->getMu(0);
        }
        FemMat.updateMat(mu);

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
            preconditionedSSORConjugateGradientSolver(FemMat, b, A, 1.5, 1e-6, 100000);
            // preconditionedIncompleteCholeskyConjugateGradientSolver(FemMat, b, A, 1e-6, 100000);
        } else {
            printf("nonlinear materials detected, starting non linear newton solver\n");
            CV r(b.val.size());  // residual

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
                FemMat.updateMu(props, mu, B, residual, i);
                FemMat.updateMat(mu);
                CV::mult(r, FemMat, A);
                CV::sub(r, b, r);
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

        auto start9 = std::chrono::high_resolution_clock::now();

        mesh.computeCurl(B, A);

        printf("%f translate and fix mesh\n", std::chrono::duration_cast<std::chrono::duration<float>>(start2 - start1).count()*1000);
        printf("%f mesh\n", std::chrono::duration_cast<std::chrono::duration<float>>(start3 - start2).count()*1000);
        printf("%f kelvin boundary conditions\n", std::chrono::duration_cast<std::chrono::duration<float>>(start4 - start3).count()*1000);
        printf("%f allocate vectors b x and coo\n", std::chrono::duration_cast<std::chrono::duration<float>>(start5 - start4).count()*1000);
        printf("%f get fem system\n", std::chrono::duration_cast<std::chrono::duration<float>>(start6 - start5).count()*1000);
        printf("%f add dirichlet boundary conditions\n", std::chrono::duration_cast<std::chrono::duration<float>>(start7 - start6).count()*1000);
        printf("%f convert to csr\n", std::chrono::duration_cast<std::chrono::duration<float>>(start8 - start7).count()*1000);
        printf("%f solve\n", std::chrono::duration_cast<std::chrono::duration<float>>(start9 - start8).count()*1000);
        printf("%f total\n", std::chrono::duration_cast<std::chrono::duration<float>>(start9 - start1).count()*1000);

        // auto integrals = mesh.computeForceIntegrals(B);
        // for (auto& i : integrals) {
        //     printf("force integral: %.17gx + %.17gy\n", i.x, i.y);
        // }
        // mesh.muplot(mu);
        // return;
        /*
        mesh.Bplot();

        // report(&out, 1, 1, 0, 0, 0, 0);

        printf("Number of points: %d\nNumber of boundary vertices: %d\n", mesh.vertices.size(), mesh.boundary_vertices.size());
        */

    }
}