#ifndef NIK_MAGNETOSTATICMESH_HPP
#define NIK_MAGNETOSTATICMESH_HPP

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <math.h>
#include <set>

#include "../../lib/triangle/triangle.h"
#include "../triangle/util.h"

#include "../constants.hpp"

#include "../utils/utils.hpp"
#include "../drawing/drawing.hpp"
#include "../mesh/mesh.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/build_coo.hpp"
#include "../constants.hpp"
#include "magnetostatic_algebra.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticMesh : Mesh<MagnetostaticProp> {
        MagnetostaticMesh() {
            // default material property
            default_prop = {0, {0, 0}, materials::air, {}, 0, {0, 0}};
        }

        ~MagnetostaticMesh() {
            
        }

        std::vector<double> getFemSystem(BuildMatCOO<MagnetostaticNonLinearExpression>& coo, CV& b);
        void addDirichletBoundaryConditions(BuildMatCOO<MagnetostaticNonLinearExpression>& coo, CV& b);
        void computeCurl(std::vector<Vector>& B, CV& A);
        void refineMeshAroundMagnets();
        std::vector<Vector> computeForceIntegrals(std::vector<Vector>& B);
    };

    std::vector<double> MagnetostaticMesh::getFemSystem(BuildMatCOO<MagnetostaticNonLinearExpression>&coo, CV &b) {
        // since the stiffness matrix is symmetric, this function only computes the upper triangular part
        auto start = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG_PRINT
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            printf("%d %d\n", i, (int)drawing.getRegionFromId(data.triangleattributelist[i]).mu);
        }
        #endif

        auto eadjelems_ids = std::vector<std::array<uint32_t, 3>>(data.numberoftriangles);
        auto eadjelems_props = std::vector<std::array<const MagnetostaticProp*, 3>>(data.numberoftriangles);
        std::vector<double> Jm(data.numberoftriangles);
        
        struct EAdjElem {
            uint8_t size = 0;
            uint32_t ids[2];
        };
        std::unordered_map<uint64_t, EAdjElem> eadjmap((data.numberofpoints - 2) * 3);

        // get the edge adjacent elements
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t v1 = data.trianglelist[i].verts[j];
                uint32_t v2 = data.trianglelist[i].verts[(j + 1) % 3];

                uint64_t e = (((uint64_t)v1 << 32) | v2) * (v1 > v2) + (((uint64_t)v2 << 32) | v1) * (v1 < v2);

                eadjmap[e].ids[eadjmap[e].size++] = i;
            }
        }

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t v1 = data.trianglelist[i].verts[j];
                uint32_t v2 = data.trianglelist[i].verts[(j + 1) % 3];

                uint64_t e = (((uint64_t)v1 << 32) | v2) * (v1 > v2) + (((uint64_t)v2 << 32) | v1) * (v1 < v2);

                eadjelems_ids[i][j] = eadjmap[e].ids[0] * (eadjmap[e].ids[0] != i) + eadjmap[e].ids[1] * (eadjmap[e].ids[1] != i);
                eadjelems_props[i][j] = drawing.getRegionPtrFromId(data.triangleattributelist[eadjelems_ids[i][j]]);
            }
        }

        // compute curl
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            // get barycenter of adjacent elements
            Point barycenters[3];
            for (uint8_t j = 0; j < 3; j++) {
                barycenters[j] = (data.pointlist[data.trianglelist[eadjelems_ids[i][j]].verts[0]] +
                                  data.pointlist[data.trianglelist[eadjelems_ids[i][j]].verts[1]] +
                                  data.pointlist[data.trianglelist[eadjelems_ids[i][j]].verts[2]]) / 3.f;
            }

            Point barycenter = (data.pointlist[data.trianglelist[i].verts[0]] +
                                data.pointlist[data.trianglelist[i].verts[1]] +
                                data.pointlist[data.trianglelist[i].verts[2]]) / 3.f;

            // get difference of magnetization vector components, adjmag - mag
            const MagnetostaticProp* mag = drawing.getRegionPtrFromId(data.triangleattributelist[i]);

            double dx1 = barycenters[0].x - barycenter.x;
            double dy1 = barycenters[0].y - barycenter.y;
            double dx2 = barycenters[1].x - barycenter.x;
            double dy2 = barycenters[1].y - barycenter.y;
            double dx3 = barycenters[2].x - barycenter.x;
            double dy3 = barycenters[2].y - barycenter.y;

            double dfx1 = eadjelems_props[i][0]->M.x - mag->M.x;
            double dfy1 = eadjelems_props[i][0]->M.y - mag->M.y;
            double dfx2 = eadjelems_props[i][1]->M.x - mag->M.x;
            double dfy2 = eadjelems_props[i][1]->M.y - mag->M.y;
            double dfx3 = eadjelems_props[i][2]->M.x - mag->M.x;
            double dfy3 = eadjelems_props[i][2]->M.y - mag->M.y;

            // fit for Mx
            double Mxa1 = dx2 - dx1;
            double Mxb1 = dy2 - dy1;
            double Mxc1 = dfx2 - dfx1;
            double Mxa2 = dx3 - dx1;
            double Mxb2 = dy3 - dy1;
            double Mxc2 = dfx3 - dfx1;

            double Mxa = Mxb1 * Mxc2 - Mxb2 * Mxc1;
            double Mxb = Mxa2 * Mxc1 - Mxa1 * Mxc2;
            double Mxc = Mxa1 * Mxb2 - Mxb1 * Mxa2;

            double Mxdx = -Mxa / Mxc;
            double Mxdy = -Mxb / Mxc;

            // fit for My
            double Mya1 = dx2 - dx1;
            double Myb1 = dy2 - dy1;
            double Myc1 = dfy2 - dfy1;
            double Mya2 = dx3 - dx1;
            double Myb2 = dy3 - dy1;
            double Myc2 = dfy3 - dfy1;

            double Mya = Myb1 * Myc2 - Myb2 * Myc1;
            double Myb = Mya2 * Myc1 - Mya1 * Myc2;
            double Myc = Mya1 * Myb2 - Myb1 * Mya2;

            double Mydx = -Mya / Myc;
            double Mydy = -Myb / Myc;

            // curl
            double curl = Mxdy - Mydx;

            Jm[i] = curl;

            /*
            Jm[i] = (dfx1*(dx1*(dx1*dy1 + dx2*dy2 + dx3*dy3) -
                    dy1*(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2))) +
                    dfx2*(dx2*(dx1*dy1 + dx2*dy2 + dx3*dy3) -
                    dy2*(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2))) +
                    dfx3*(dx3*(dx1*dy1 + dx2*dy2 + dx3*dy3) -
                    dy3*(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2))) +
                    dfy1*(dx1*(pow(dy1, 2) + pow(dy2, 2) + pow(dy3, 2)) -
                    dy1*(dx1*dy1 + dx2*dy2 + dx3*dy3)) +
                    dfy2*(dx2*(pow(dy1, 2) + pow(dy2, 2) + pow(dy3, 2)) -
                    dy2*(dx1*dy1 + dx2*dy2 + dx3*dy3)) +
                    dfy3*(dx3*(pow(dy1, 2) + pow(dy2, 2) + pow(dy3, 2)) - 
                    dy3*(dx1*dy1 + dx2*dy2 + dx3*dy3)))/(pow(dx1, 2)*pow(dy2, 2) + 
                    pow(dx1, 2)*pow(dy3, 2) - 2*dx1*dx2*dy1*dy2 - 
                    2*dx1*dx3*dy1*dy3 + pow(dx2, 2)*pow(dy1, 2) + 
                    pow(dx2, 2)*pow(dy3, 2) - 2*dx2*dx3*dy2*dy3 + 
                    pow(dx3, 2)*pow(dy1, 2) + pow(dx3, 2)*pow(dy2, 2));
            */
            
            // printf("Jm[%d] = %.17g\n", i, Jm[i]);
        }

        /*
        auto adjelems_ids = new uint32_t[data.numberofpoints][18];
        auto adjelems_props = new const MagnetostaticProp*[data.numberofpoints][18];
        auto adjelems_count = new uint8_t[data.numberofpoints]();  // initialize to 0
        */

        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(data.numberofpoints);
        auto adjelems_props = std::vector<std::array<const MagnetostaticProp*, 18>>(data.numberofpoints);
        auto adjelems_count = std::vector<uint8_t>(data.numberofpoints, 0);

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = data.trianglelist[i].verts[j];
                adjelems_ids[myid][adjelems_count[myid]] = i;
                adjelems_props[myid][adjelems_count[myid]++] = drawing.getRegionPtrFromId(data.triangleattributelist[i]);
            }
        }
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            for (uint8_t j = 0; j < adjelems_count[i]; j++) {
                uint32_t v1, v2, v3;
                v1 = i;
                Elem myelem = data.trianglelist[adjelems_ids[i][j]];
                if (i == data.trianglelist[adjelems_ids[i][j]].verts[0]) {
                    v2 = myelem.verts[1];
                    v3 = myelem.verts[2];
                } else if (i == myelem.verts[1]) {
                    v2 = myelem.verts[2];
                    v3 = myelem.verts[0];
                } else if (i == myelem.verts[2]) {
                    v2 = myelem.verts[0];
                    v3 = myelem.verts[1];
                } else {
                    nexit("error: vertex not found in element");
                }

                double oriented_area = Point::doubleOrientedArea(data.pointlist[v1], data.pointlist[v2], data.pointlist[v3]);
                
                if (oriented_area < 0) {
                    std::swap(v2, v3);
                }

                double area = Point::doubleOrientedArea(data.pointlist[v1], data.pointlist[v2], data.pointlist[v3]);

                // elements are only added if they are in the upper triangle because the matrix is symmetric and this saves half the memory
                double b1 = (data.pointlist[v2].y - data.pointlist[v3].y) / area;
                double c1 = (data.pointlist[v3].x - data.pointlist[v2].x) / area;
                // coo.add_elem(i, v1, (area * (b1 * b1 + c1 * c1)) / (2 * adjelems_props[i][j].mu));
                // check if key exists
                if (v1 <= i) {
                    uint64_t key = BuildMatCOO<int>::get_key(i, v1);
                    coo.elems.emplace(key, MagnetostaticNonLinearExpression());
                    coo.elems[key].terms.push_back(
                        {
                            (area * (b1 * b1 + c1 * c1) * 0.5),
                            adjelems_ids[i][j],
                            false
                        }
                    );
                }
                if (v2 <= i) {
                    double b2 = (data.pointlist[v3].y - data.pointlist[v1].y) / area;
                    double c2 = (data.pointlist[v1].x - data.pointlist[v3].x) / area;
                    // coo.add_elem(i, v2, (area * (b2 * b1 + c2 * c1)) / (2 * adjelems_props[i][j].mu));
                    uint64_t key = BuildMatCOO<int>::get_key(i, v2);
                    coo.elems.emplace(key, MagnetostaticNonLinearExpression());
                    coo.elems[key].terms.push_back(
                        {
                            (area * (b2 * b1 + c2 * c1) * 0.5),
                            adjelems_ids[i][j],
                            false
                        }
                    );
                }
                if (v3 <= i) {
                    double b3 = (data.pointlist[v1].y - data.pointlist[v2].y) / area;
                    double c3 = (data.pointlist[v2].x - data.pointlist[v1].x) / area;
                    // coo.add_elem(i, v3, (area * (b3 * b1 + c3 * c1)) / (2 * adjelems_props[i][j].mu));
                    uint64_t key = BuildMatCOO<int>::get_key(i, v3);
                    coo.elems.emplace(key, MagnetostaticNonLinearExpression());
                    coo.elems[key].terms.push_back(
                        {
                            (area * (b3 * b1 + c3 * c1) * 0.5),
                            adjelems_ids[i][j],
                            false
                        }
                    );
                }

                // set the b vector
                b.add_elem(i, (area * (adjelems_props[i][j]->J + Jm[adjelems_ids[i][j]])) / 6);
            }
        }

        // iterate over upper triangular matrix and copy to lower triangular matrix
        // for (uint32_t i = 0; i < coo.m; i++) {
        //     for (uint32_t j = i; j < coo.n; j++) {
        //         if (coo.get_elem(i, j) != 0) {
        //             coo.add_elem(j, i, coo.get_elem(i, j));
        //         }
        //     }
        // }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "FEM matrix construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
        return Jm;
    }

    void MagnetostaticMesh::addDirichletBoundaryConditions(BuildMatCOO<MagnetostaticNonLinearExpression> &coo, CV &b) {
        // find three furthest points from the center
        uint32_t p1, p2, p3;
        double d1, d2, d3;
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            double dist = Point::distance(data.pointlist[i], Point(0, 0));
            if (dist > d1) {
                d3 = d2;
                p3 = p2;
                d2 = d1;
                p2 = p1;
                d1 = dist;
                p1 = i;
            } else if (dist > d2) {
                d3 = d2;
                p3 = p2;
                d2 = dist;
                p2 = i;
            } else if (dist > d3) {
                d3 = dist;
                p3 = i;
            }
        }

        for (auto elem : coo.elems) {
            uint32_t m = elem.first >> 32;
            uint32_t n = elem.first & 0xFFFFFFFF;
            if (m == p1 || m == p2 || m == p3) {
                coo.elems[BuildMatCOO<int>::get_key(m, n)].terms.clear();
                coo.elems[BuildMatCOO<int>::get_key(m, n)].terms.push_back({0, 0, true});
            }
            if (n == p1 || n == p2 || n == p3) {
                coo.elems[BuildMatCOO<int>::get_key(m, n)].terms.clear();
                coo.elems[BuildMatCOO<int>::get_key(m, n)].terms.push_back({0, 0, true});
            }
        }

        coo.elems[BuildMatCOO<int>::get_key(p1, p1)].terms.clear();
        coo.elems[BuildMatCOO<int>::get_key(p1, p1)].terms.push_back({ 1, 0, true});
        coo.elems[BuildMatCOO<int>::get_key(p2, p2)].terms.clear();
        coo.elems[BuildMatCOO<int>::get_key(p2, p2)].terms.push_back({ 1, 0, true});
        coo.elems[BuildMatCOO<int>::get_key(p3, p3)].terms.clear();
        coo.elems[BuildMatCOO<int>::get_key(p3, p3)].terms.push_back({ 1, 0, true});
        // since we are setting the potential to zero at the three points, we do not need to subtract anything from the b vector
    }

    void MagnetostaticMesh::computeCurl(std::vector<Vector>& B, CV &A) {
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            Elem myelem = data.trianglelist[i];
            double x1 = data.pointlist[myelem.verts[0]].x;
            double y1 = data.pointlist[myelem.verts[0]].y;
            double z1 = A[myelem.verts[0]];
            double x2 = data.pointlist[myelem.verts[1]].x;
            double y2 = data.pointlist[myelem.verts[1]].y;
            double z2 = A[myelem.verts[1]];
            double x3 = data.pointlist[myelem.verts[2]].x;
            double y3 = data.pointlist[myelem.verts[2]].y;
            double z3 = A[myelem.verts[2]];

            // fit z = a + bx + cy to the three points
            double a1 = x2 - x1;
            double b1 = y2 - y1;
            double c1 = z2 - z1;
            double a2 = x3 - x1;
            double b2 = y3 - y1;
            double c2 = z3 - z1;

            double a = b1 * c2 - b2 * c1;
            double b = a2 * c1 - a1 * c2;
            double c = a1 * b2 - b1 * a2;

            double dx = -a / c;
            double dy = -b / c;
            
            // printf("dx = %.17g dy = %.17g for elem (%.1f, %.1f, %.17g), (%.1f, %.1f, %.17g), (%.1f, %.1f, %.17g)\n", dx, dy, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            B[i] = Vector(-dy, dx);
        }
    }

    void MagnetostaticMesh::refineMeshAroundMagnets() {

        // add magnets refining points
        // first we have to find the polygons that are magnets
        // if the polygon contains a region that is a magnet and it does not contain any other polygon that contains the same region then it is a magnet
        std::map<uint32_t, std::vector<Polygon>> polygons_that_contain_magnet_region;
        for (auto region :drawing.regions) {
            if (drawing.region_map[region.second].M != Vector(0, 0)) {
                for (auto polygon : drawing.polygons) {
                    if (polygon.contains(region.first)) {
                        polygons_that_contain_magnet_region[region.second].push_back(polygon);
                    }
                }
            }
        }
        std::vector<Polygon> magnets;
        for (auto region_poly_array : polygons_that_contain_magnet_region) {
            if (region_poly_array.second.size() == 1) {
                // this is a magnet
                for (auto point : region_poly_array.second[0].points) {
                    magnets.push_back(region_poly_array.second[0]);
                }
            } else {
                // this array contains a magnet and other polygons that contain the polygon that is a magnet
                // we need to find the polygon that does not contain any other polygon that contains the same region
                for (auto polygon : region_poly_array.second) {
                    bool is_magnet = true;
                    for (auto other_polygon : region_poly_array.second) {
                        if (polygon != other_polygon && polygon.contains(other_polygon)) {
                            is_magnet = false;
                            break;
                        }
                    }
                    if (is_magnet) {
                        magnets.push_back(polygon);
                    }
                }
            }
        }
        printf("Found %lu magnets\n", magnets.size());
        for (auto polygon : magnets) {
            for (uint32_t i = 0; i < polygon.points.size(); i++) {
                Point p1 = polygon.points[i];
                Point p2 = polygon.points[(i + 1) % polygon.points.size()];
                
                double multiplier = 10;
                uint32_t n_points = Point::distance(p1, p2) / (drawing.epsilon / multiplier);
                // multiply by 1.1 to make sure that the points are far enough to not trigger conformity checks
                Vector normal = Vector::normal(p1, p2).versor() * (drawing.epsilon / multiplier) * 1.1;

                for (uint32_t j = 0; j < n_points; j++) {
                    Point mypoints[2];
                    mypoints[0] = Point::lerp(p1, p2, (double)j / n_points) + normal;
                    mypoints[1] = Point::lerp(p1, p2, (double)j / n_points) - normal;

                    bool mypoints_ok[2] = {
                        polygon.contains(mypoints[0]), 
                        polygon.contains(mypoints[1])
                    };

                    // check points at least drawing.epsilon away from each edge of each polygon
                    for (uint8_t k = 0; k < 2; k++) {
                        // skip check for failed points
                        if (mypoints_ok[k]) {
                            // for every polygon
                            for (auto other_polygon : drawing.polygons) {
                                // for every edge of the other polygon
                                for (uint32_t l = 0; l < other_polygon.points.size(); l++) {
                                    // check if the point is at least drawing.epsilon away from the edge
                                    Point p3 = other_polygon.points[l];
                                    Point p4 = other_polygon.points[(l + 1) % other_polygon.points.size()];
                                    if (Point::distance(mypoints[k], p3) < drawing.epsilon / multiplier ||
                                        Point::distance(mypoints[k], p4) < drawing.epsilon / multiplier) {
                                        mypoints_ok[k] = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    // check points at least drawing.epsilon away from each point in drawing.points
                    for (uint8_t k = 0; k < 2; k++) {
                        // skip check for failed points
                        if (mypoints_ok[k]) {
                            for (auto point : drawing.points) {
                                if (Point::distance(mypoints[k], point) < drawing.epsilon / multiplier) {
                                    mypoints_ok[k] = false;
                                    break;
                                }
                            }
                        }
                    }

                    // add points to drawing.points
                    for (uint8_t k = 0; k < 2; k++) {
                        if (mypoints_ok[k]) {
                            drawing.points.push_back(mypoints[k]);
                        }
                    }
                }
            }
        }
    }

    std::vector<Vector> MagnetostaticMesh::computeForceIntegrals(std::vector<Vector>& B) {
        std::vector<Vector> force_integrals(drawing.polygons.size());
        // make epsilon a little smaller than the smallest distance
        for (auto polygon : drawing.polygons) {
            force_integrals.push_back(Vector(0, 0));
            for (uint32_t i = 0; i < polygon.points.size(); i++) {
                Point p1 = polygon.points[i];
                Point p2 = polygon.points[(i + 1) % polygon.points.size()];
                printf("p1 = (%.17g, %.17g) p2 = (%.17g, %.17g)\n", p1.x, p1.y, p2.x, p2.y);
                // evaluate integral over line segment
                Vector normal = Vector(p2.y - p1.y, p1.x - p2.x) / ((Vector)(p2 - p1)).magnitude();
                std::vector<Point> integration_segment_mesh_intersection_points;
                // find vertices that are on the line segment
                for (uint32_t j = 0; j < data.numberofpoints; j++) {
                    Point p = data.pointlist[j];
                    double dist = Segment::pointSegmentDistance(p, p1, p2);
                    if (dist < epsilon) {
                        integration_segment_mesh_intersection_points.push_back(p);
                    }
                }
                printf("size %d\n", integration_segment_mesh_intersection_points.size());
                // sort vertices by distance from p1
                std::sort(integration_segment_mesh_intersection_points.begin(), integration_segment_mesh_intersection_points.end(), [&p1](Point a, Point b) {
                    return Point::distance(a, p1) < Point::distance(b, p1);
                });
                // find adj vertices to vertices that are on the line segment
                std::vector<std::vector<Vector>> adj_values(integration_segment_mesh_intersection_points.size(), std::vector<Vector>()); 
                for (uint32_t e = 0; e < data.numberoftriangles; e++) {
                    Elem myelem = data.trianglelist[e];
                    for (uint32_t j = 0; j < 3; j++) {
                        Point p = data.pointlist[myelem.verts[j]];
                        for (uint32_t k = 0; k < integration_segment_mesh_intersection_points.size(); k++) {
                            if (Point::distance(p, integration_segment_mesh_intersection_points[k]) < epsilon) {
                                adj_values[k].push_back(B[e]);
                            }
                        }
                    }
                }
                // this could be a least squares fit but it's probably not necessary for now
                std::vector<Vector> adj_values_avg(integration_segment_mesh_intersection_points.size(), {0, 0});
                for (uint32_t i = 0; i < adj_values.size(); i++) {
                    for (uint32_t j = 0; j < adj_values[i].size(); j++) {
                        adj_values_avg[i] += adj_values[i][j];
                    }
                    adj_values_avg[i] /= adj_values[i].size();
                }
                // trapezoidal integration
                for (uint32_t i = 0; i < integration_segment_mesh_intersection_points.size() - 1; i++) {
                    Point p1 = integration_segment_mesh_intersection_points[i];
                    Point p2 = integration_segment_mesh_intersection_points[i + 1];
                    Vector v1 = adj_values_avg[i];
                    Vector v2 = adj_values_avg[i + 1];
                    MaxwellStressTensor stress_tensor1;
                    // 1 / mu_0 * (BiBj - 1/2 * B*B * delta_ij)
                    stress_tensor1[0][0] = (1. / MU_0) * (v1.x * v1.x - 0.5 * v1.magnitude() * v1.magnitude());
                    stress_tensor1[0][1] = (1. / MU_0) * (v1.x * v1.y);
                    stress_tensor1[1][0] = (1. / MU_0) * (v1.y * v1.x);
                    stress_tensor1[1][1] = (1. / MU_0) * (v1.y * v1.y - 0.5 * v1.magnitude() * v1.magnitude());

                    MaxwellStressTensor stress_tensor2;
                    stress_tensor2[0][0] = (1. / MU_0) * (v2.x * v2.x - 0.5 * v2.magnitude() * v2.magnitude());
                    stress_tensor2[0][1] = (1. / MU_0) * (v2.x * v2.y);
                    stress_tensor2[1][0] = (1. / MU_0) * (v2.y * v2.x);
                    stress_tensor2[1][1] = (1. / MU_0) * (v2.y * v2.y - 0.5 * v2.magnitude() * v2.magnitude());

                    Vector force1(stress_tensor1[0][0] * normal.x + stress_tensor1[0][1] * normal.y, stress_tensor1[1][0] * normal.x + stress_tensor1[1][1] * normal.y);
                    Vector force2(stress_tensor2[0][0] * normal.x + stress_tensor2[0][1] * normal.y, stress_tensor2[1][0] * normal.x + stress_tensor2[1][1] * normal.y);

                    Vector force = (force1 + force2) * 0.5 * Point::distance(p1, p2);
                    printf("point1 = (%.17g, %.17g) point2 = (%.17g, %.17g) force = (%.17g, %.17g)\n", p1.x, p1.y, p2.x, p2.y, force.x, force.y);
                    force_integrals.back() += force;
                }
            }
        }
        return force_integrals;
    }
}

#endif