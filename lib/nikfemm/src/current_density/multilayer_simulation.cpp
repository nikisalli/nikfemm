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

#include "../../lib/triangle/triangle.h"

#include "../constants.hpp"
#include "multilayer_simulation.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/vector.hpp"
#include "../algebra/build_coo.hpp"
#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/solvers.hpp"

namespace nikfemm {
    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers, std::vector<double> depths, std::vector<double> max_triangle_areas) {
        if (depths.size() != num_layers) {
            throw std::invalid_argument("depths.size() != num_layers");
        }
        if (max_triangle_areas.size() != num_layers) {
            throw std::invalid_argument("max_triangle_areas.size() != num_layers");
        }

        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].max_triangle_area = max_triangle_areas[i];
            meshes[i].depth = depths[i];
        }
    }

    MultiLayerCurrentDensitySimulation::MultiLayerCurrentDensitySimulation(uint32_t num_layers) {
        for (uint32_t i = 0; i < num_layers; i++) {
            meshes.push_back(CurrentDensityMesh());
            meshes[i].max_triangle_area = 1;
            meshes[i].depth = 1;
        }
    }

    MultiLayerCurrentDensitySimulation::~MultiLayerCurrentDensitySimulation() {

    }

    static inline double geomAngle(Vector a, Vector b, Vector c) {
        // safe angle calculation
        double ax = a.x - b.x;
        double ay = a.y - b.y;
        double bx = c.x - b.x;
        double by = c.y - b.y;
        double dot = ax * bx + ay * by;
        double det = ax * by - ay * bx;
        double angle = fabs(atan2(det, dot));
        return angle;
    }

    CurrentDensitySystem MultiLayerCurrentDensitySimulation::generateSystem(bool refine) {
        /*
        if (refine) {
            for (auto& mesh : meshes) {
                mesh.drawing.addRefiningPoints();
            }
        }

        // method of the wormholes (the vias are wormholes that connect two different 2d layers)
        // we will create tunnels in our meshes that simulate the vias
        for (auto interconnection : interconnections) {
            meshes[interconnection.layer1_id].drawing.drawCircle(interconnection.p1, interconnection.diameter / 2, 16);
            meshes[interconnection.layer2_id].drawing.drawCircle(interconnection.p2, interconnection.diameter / 2, 16);
        }

        // mesh the layers
        for (uint64_t i = 0; i < meshes.size(); i++) {
            meshes[i].mesh();
            meshes[i].computeEpsilon();
            nloginfo("Meshed layer %d with %d nodes and %d elements", i, meshes[i].data.numberofpoints, meshes[i].data.numberoftriangles);
        }

        // now we need to manually merge the meshes into one big mesh and connect the circles to one another by adding triangles
        auto mymesh = CurrentDensityMesh();
        // add the points from all the meshes
        uint64_t total_points = 0;
        for (auto& mesh : meshes) {
            total_points += mesh.data.numberofpoints;
        }
        // reallocate the pointlist
        mymesh.data.pointlist = (Vector*) realloc(mymesh.data.pointlist, total_points * sizeof(Vector));
        // now add the points from all the meshes
        uint64_t offset = 0;
        for (auto& mesh : meshes) {
            for (uint64_t i = 0; i < mesh.data.numberofpoints; i++) {
                mymesh.data.pointlist[i + offset] = mesh.data.pointlist[i];
            }
            offset += mesh.data.numberofpoints;
        }

        // add the triangles from all the meshes
        uint64_t total_triangles = 0;
        for (auto& mesh : meshes) {
            total_triangles += mesh.data.numberoftriangles;
        }
        // reallocate the trianglelist
        mymesh.data.trianglelist = (Elem*) realloc(mymesh.data.trianglelist, total_triangles * sizeof(Elem));
        // now add the triangles from all the meshes
        offset = 0;
        for (auto& mesh : meshes) {
            for (uint64_t i = 0; i < mesh.data.numberoftriangles; i++) {
                mymesh.data.trianglelist[i + offset] = mesh.data.trianglelist[i];
            }
            offset += mesh.data.numberoftriangles;
        }

        // add the triangle attributes from all the meshes
        uint64_t total_triangle_attributes = 0;
        for (auto& mesh : meshes) {
            total_triangle_attributes += mesh.data.numberoftriangleattributes;
        }

        // build some useful offsets
        std::vector<uint64_t> region_ids_offsets;
        uint64_t region_offset = 0;
        for (auto& mesh : meshes) {
            region_ids_offsets.push_back(region_offset);
            region_offset += mesh.drawing.region_map.size();
        }

        // reallocate the triangleattributelist
        mymesh.data.triangleattributelist = (TRI_REAL*) realloc(mymesh.data.triangleattributelist, total_triangle_attributes * sizeof(TRI_REAL));
        // now add the triangle attributes from all the meshes
        offset = 0;
        for (uint64_t i = 0; i < meshes.size(); i++) {
            for (uint64_t j = 0; j < meshes[i].data.numberoftriangleattributes; j++) {
                mymesh.data.triangleattributelist[j + offset] = meshes[i].data.triangleattributelist[j] + region_ids_offsets[i];
            }
            offset += meshes[i].data.numberoftriangleattributes;
        }

        // we need an updated region map merging all the region maps from the meshes
        mymesh.drawing.region_map.clear();
        for (uint64_t i = 0; i < meshes.size(); i++) {
            for (auto& region : meshes[i].drawing.region_map) {
                mymesh.drawing.region_map.push_back(region);
                // we may want to deduplicate this map but why even bother?? it's not gonna give us any performance boost
            }
        }

        // we only need to add points, triangles and triangle attributes because all the other
        // fields are f**ked (censored because copilot won't autocomplete otherwise) since
        // our mesh is some kind of 2d manifold full of holes

        // now we need to actually connect the holes to one another

        // build an offset list so we can easily navigate the list of points of our manifold
        std::vector<uint64_t> offsets;
        offset = 0;
        for (auto& mesh : meshes) {
            offsets.push_back(offset);
            offset += mesh.data.numberofpoints;
        }

        // let's build a list of useful structures
        struct hole {
            CurrentDensityInterconnection interconnection;
            std::vector<uint64_t> layer1_points;
            std::vector<uint64_t> layer2_points;
        };

        std::vector<hole> holes;
        
        // for every interconnection we need to find the points that are inside the circles using the precomputed epsilon
        for (auto& interconnection : interconnections) {
            // find the points that are inside the circles
            std::vector<uint64_t> layer1_points;
            std::vector<uint64_t> layer2_points;

            auto& mesh1 = meshes[interconnection.layer1_id];
            auto& mesh2 = meshes[interconnection.layer2_id];

            for (uint64_t i = 0; i < mesh1.data.numberofpoints; i++) {
                if (fabs((mesh1.data.pointlist[i] - interconnection.p1).norm() - interconnection.diameter / 2) < mesh1.epsilon) {
                    layer1_points.push_back(i);
                }
            }

            for (uint64_t i = 0; i < mesh2.data.numberofpoints; i++) {
                if (fabs((mesh2.data.pointlist[i] - interconnection.p2).norm() - interconnection.diameter / 2) < mesh2.epsilon) {
                    layer2_points.push_back(i);
                }
            }

            // print some info
            nloginfo("found %d points in layer %d that are inside the circle", layer1_points.size(), interconnection.layer1_id);
            nloginfo("found %d points in layer %d that are inside the circle", layer2_points.size(), interconnection.layer2_id);

            holes.push_back({interconnection, layer1_points, layer2_points});
        }

        // now we need to connect the holes to one another by adding triangles
        // we will do this by having a list of equivalent points and connecting them in couples
        // we will also need to add the triangles to the triangleattributelist

        std::vector<std::pair<Elem, uint64_t>> triangles_to_add;
        for (auto& hole : holes) {
            // find the equivalent points by comparing the distances using the smallest epsilon between the two layers
            auto& mesh1 = meshes[hole.interconnection.layer1_id];
            auto& mesh2 = meshes[hole.interconnection.layer2_id];

            // find the smallest epsilon
            double epsilon = mesh1.epsilon;
            if (mesh2.epsilon < epsilon) {
                epsilon = mesh2.epsilon;
            }

            nloginfo("epsilon: %f", epsilon);

            // find the equivalent points
            // first point is always on the first layer and the second point is always on the second layer
            std::vector<std::pair<uint64_t, uint64_t>> equivalent_points;
            for (auto& point1 : hole.layer1_points) {
                for (auto& point2 : hole.layer2_points) {
                    if (((mesh1.data.pointlist[point1] - hole.interconnection.p1) - (mesh2.data.pointlist[point2] - hole.interconnection.p2)).norm() < epsilon) {
                        equivalent_points.push_back({point1, point2});
                    }
                }
            }

            // print some info
            nloginfo("found %d equivalent points", equivalent_points.size());
            for (auto& point : equivalent_points) {
                nloginfo("point: %d, %d", point.first, point.second);
            }

            // the hole is essentially a prism, we now need to find all the lateral faces of the prism
            // a lateral face is essentially a rectangle formed by two couples of equivalent points
            // we will use a set to avoid duplicates
            std::set<std::pair<uint64_t, uint64_t>> lateral_faces;  // the indices here are the indices of the equivalent points couples
            std::set<std::pair<uint64_t, uint64_t>> adj_points1; // adjacent points to the first point (layer 1)
            std::set<std::pair<uint64_t, uint64_t>> adj_points2; // adjacent points to the second point (layer 2)
            // find the adjacent points
            for (auto& point : equivalent_points) {
                for (uint64_t j = 0; j < mesh1.data.numberoftriangles; j++) {
                    auto& triangle = mesh1.data.trianglelist[j];
                    for (uint64_t i = 0; i < 3; i++) {
                        if (triangle[i] == point.first) {
                            // only one of the two adjacent points will be in the equivalent points list
                            // find which one it is
                            bool p1_found = false;
                            bool p2_found = false;
                            for (auto& p : equivalent_points) {
                                if (p.first == triangle[(i + 1) % 3]) {
                                    p1_found = true;
                                }
                                if (p.first == triangle[(i + 2) % 3]) {
                                    p2_found = true;
                                }
                            }
                            if (p1_found && p2_found) {
                                // should never happen
                                nexit("error: found two adjacent points in the equivalent points list");
                            } else if (p1_found) {
                                // insert but sort the pair
                                if (point.first < triangle[(i + 1) % 3]) {
                                    adj_points1.insert(std::pair<uint64_t, uint64_t>(point.first, triangle[(i + 1) % 3]));
                                } else {
                                    adj_points1.insert(std::pair<uint64_t, uint64_t>(triangle[(i + 1) % 3], point.first));
                                }
                            } else if (p2_found) {
                                // insert but sort the pair
                                if (point.first < triangle[(i + 2) % 3]) {
                                    adj_points1.insert(std::pair<uint64_t, uint64_t>(point.first, triangle[(i + 2) % 3]));
                                } else {
                                    adj_points1.insert(std::pair<uint64_t, uint64_t>(triangle[(i + 2) % 3], point.first));
                                }
                            } else {
                                // only one point in this triangle, skip
                            }
                        }
                    }
                }

                for (uint64_t j = 0; j < mesh2.data.numberoftriangles; j++) {
                    auto& triangle = mesh2.data.trianglelist[j];
                    for (uint64_t i = 0; i < 3; i++) {
                        if (triangle[i] == point.second) {
                            // only one of the two adjacent points will be in the equivalent points list
                            // find which one it is
                            bool p1_found = false;
                            bool p2_found = false;
                            for (auto& p : equivalent_points) {
                                if (p.second == triangle[(i + 1) % 3]) {
                                    p1_found = true;
                                }
                                if (p.second == triangle[(i + 2) % 3]) {
                                    p2_found = true;
                                }
                            }
                            if (p1_found && p2_found) {
                                // should never happen
                                nexit("error: found two adjacent points in the equivalent points list");
                            } else if (p1_found) {
                                // insert but sort the pair
                                if (point.second < triangle[(i + 1) % 3]) {
                                    adj_points2.insert(std::pair<uint64_t, uint64_t>(point.second, triangle[(i + 1) % 3]));
                                } else {
                                    adj_points2.insert(std::pair<uint64_t, uint64_t>(triangle[(i + 1) % 3], point.second));
                                }
                            } else if (p2_found) {
                                // insert but sort the pair
                                if (point.second < triangle[(i + 2) % 3]) {
                                    adj_points2.insert(std::pair<uint64_t, uint64_t>(point.second, triangle[(i + 2) % 3]));
                                } else {
                                    adj_points2.insert(std::pair<uint64_t, uint64_t>(triangle[(i + 2) % 3], point.second));
                                }
                            } else {
                                // only one point in this triangle, skip
                            }
                        }
                    }
                }
            }

            // print some info
            nloginfo("found %d adjacent points to the first point", adj_points1.size());
            for (auto& point : adj_points1) {
                nloginfo("point: %d, %d", point.first, point.second);
            }
            nloginfo("found %d adjacent points to the second point", adj_points2.size());
            for (auto& point : adj_points2) {
                nloginfo("point: %d, %d", point.first, point.second);
            }
            nloginfo("found %d lateral faces", lateral_faces.size());
            for (auto& face : lateral_faces) {
                nloginfo("face: %d, %d; %d, %d", equivalent_points[face.first].first, equivalent_points[face.first].second, equivalent_points[face.second].first, equivalent_points[face.second].second);
            }

            // now we need to add the triangles to trianglelist and triangleattributelist
            // since we have the faces we just need to choose which opposite points will make the diagonal
            // this choice is absolutely arbitrary
            for (auto& face : lateral_faces) {
                triangles_to_add.push_back({{
                    equivalent_points[face.first].first + offsets[hole.interconnection.layer1_id],
                    equivalent_points[face.first].second + offsets[hole.interconnection.layer2_id],
                    equivalent_points[face.second].second + offsets[hole.interconnection.layer2_id]
                }, mymesh.drawing.getRegionId({hole.interconnection.sigma})});
                triangles_to_add.push_back({{
                    equivalent_points[face.first].first + offsets[hole.interconnection.layer1_id],
                    equivalent_points[face.second].first + offsets[hole.interconnection.layer1_id],
                    equivalent_points[face.second].second + offsets[hole.interconnection.layer2_id]
                }, mymesh.drawing.getRegionId({hole.interconnection.sigma})});
            }
        }

        // print some info
        nloginfo("found %d triangles to add", triangles_to_add.size());
        for (auto& triangle : triangles_to_add) {
            nloginfo("triangle: %d, %d, %d", triangle.first[0], triangle.first[1], triangle.first[2]);
        }

        // now we need to add the triangles to the mesh
        // we will use the trianglelist and triangleattributelist fields
        // first we need to allocate the memory
        mymesh.data.trianglelist = (Elem*) realloc(mymesh.data.trianglelist, (mymesh.data.numberoftriangles + triangles_to_add.size()) * sizeof(Elem));
        mymesh.data.triangleattributelist = (TRI_REAL*) realloc(mymesh.data.triangleattributelist, (mymesh.data.numberoftriangleattributes + triangles_to_add.size()) * sizeof(TRI_REAL));

        // now we need to add the triangles
        for (uint64_t i = 0; i < triangles_to_add.size(); i++) {
            mymesh.data.trianglelist[mymesh.data.numberoftriangles + i] = triangles_to_add[i].first;
            mymesh.data.triangleattributelist[mymesh.data.numberoftriangleattributes + i] = triangles_to_add[i].second;
        }
        mymesh.data.numberoftriangles += triangles_to_add.size();
        mymesh.data.numberoftriangleattributes += triangles_to_add.size();

        nexit("not implemented");
        */

        if (refine) {
            for (auto& mesh : meshes) {
                mesh.drawing.addRefiningPoints();
            }
        }

        // add points where the interconnections are so that we can have the precise point in the mesh for it
        // insteaad of relying on the nearest point
        for (auto& interconnection : interconnections) {
            // the drawing algorithm will check internally if the point is already in the mesh
            // I could check if the point is inside the simulation domain but it will error out anyways if it is not
            // so I will just let it be
            meshes[interconnection.layer1_id].drawing.drawPoint(interconnection.p1);
            meshes[interconnection.layer2_id].drawing.drawPoint(interconnection.p2);
        }

        for (uint64_t i = 0; i < meshes.size(); i++) {
            meshes[i].mesh();
            meshes[i].computeEpsilon();
            nloginfo("Meshed layer %d with %d nodes and %d elements", i, meshes[i].data.numberofpoints, meshes[i].data.numberoftriangles);
        }

        std::vector<CurrentDensitySystem> systems;

        for (uint64_t i = 0; i < meshes.size(); i++) {
            systems.push_back(meshes[i].getFemSystem());
        }

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& system : systems) {
            offsets.push_back(offset);
            offset += system.b.val.size();
        }

        // now we need to merge the systems together using the interconnections
        struct interconnection_indices {
            uint64_t mesh1_id;
            uint64_t mesh2_id;
            uint64_t node1_id;
            uint64_t node2_id;
            double R;
        };
        std::vector<interconnection_indices> interconnection_indices_list;

        for (auto& interconnection : interconnections) {
            // find the nearest node in the mesh
            uint64_t node1_id = 0;
            double node1_dist = INFINITY;
            uint64_t node2_id = 0;
            double node2_dist = INFINITY;

            for (uint64_t i = 0; i < meshes[interconnection.layer1_id].data.numberofpoints; i++) {
                double dist = 0;

                dist = (meshes[interconnection.layer1_id].data.pointlist[i] - interconnection.p1).norm();
                if (dist < node1_dist) {
                    node1_id = i;
                    node1_dist = dist;
                }
            }

            for (uint64_t i = 0; i < meshes[interconnection.layer2_id].data.numberofpoints; i++) {
                double dist = 0;

                dist = (meshes[interconnection.layer2_id].data.pointlist[i] - interconnection.p2).norm();
                if (dist < node2_dist) {
                    node2_id = i;
                    node2_dist = dist;
                }
            }

            nloginfo("Nearest node to p1 is %d at x: %f, y: %f", node1_id, meshes[interconnection.layer1_id].data.pointlist[node1_id].x, meshes[interconnection.layer1_id].data.pointlist[node1_id].y);
            nloginfo("Nearest node to p2 is %d at x: %f, y: %f", node2_id, meshes[interconnection.layer2_id].data.pointlist[node2_id].x, meshes[interconnection.layer2_id].data.pointlist[node2_id].y);

            interconnection_indices_list.push_back({interconnection.layer1_id, interconnection.layer2_id, node1_id, node2_id, interconnection.R});
        }

        // preallocate the size of the merged system by summing the sizes of the systems
        uint64_t merged_system_nnz = 0;
        uint64_t merged_system_size = 0;
        for (auto& system : systems) {
            merged_system_nnz += system.A.elems.size();
            merged_system_size += system.b.val.size();
        }

        CurrentDensitySystem merged_system = {
            BuildMatCOO<double>(merged_system_size),
            CV(merged_system_size)
        };

        nloginfo("merged_system_nnz: %d", merged_system_nnz);
        nloginfo("merged_system_size: %d", merged_system_size);
        merged_system.A.elems.reserve(merged_system_nnz);

        // merge the systems together like a block matrix
        offset = 0; // reusing offset variable

        for (uint64_t i = 0; i < meshes.size(); i++) {
            CurrentDensitySystem& system = systems[i];
            CurrentDensityMesh& mesh = meshes[i];

            uint64_t system_size = mesh.data.numberofpoints;

            // add the system to the merged system
            // A
            for (auto& elem : system.A.elems) {
                uint64_t m = elem.first >> 32;
                uint64_t n = elem.first & 0xFFFFFFFF;

                merged_system.A(m + offset, n + offset) = elem.second;
            }

            // b
            for (uint64_t j = 0; j < system_size; j++) {
                merged_system.b.val[j + offset] = system.b.val[j];
            }

            offset += system_size;
        }

        nloginfo("b.size(): %d", merged_system.b.val.size());

        // add the interconnections to the merged system
        // method 1 (couldn't bother to get the math right)
        for (auto& interconnection : interconnection_indices_list) {
            // add the interconnection to the merged system
            uint64_t m1 = interconnection.node1_id + offsets[interconnection.mesh1_id];
            uint64_t m2 = interconnection.node2_id + offsets[interconnection.mesh2_id];
            double R = interconnection.R;

            // now for every row that has a nonzero element in column m1 we should:
            // 1. for every nonzero element in that row, multiply it by R and add + 1 if the column index is m1
            // 2. add a new element to the merged system with the value -1 in column m2 and the same row index

            // then for every row that has a nonzero element in column m2 we should:
            // 1. for every nonzero element in that row, multiply it by R and add + 1 if the column index is m2
            // 2. add a new element to the merged system with the value -1 in column m1 and the same row index

            // 1.
            merged_system.A(m1, m1) += 2/R;
            merged_system.A(m2, m2) += 2/R;

            // 2.
            if (m1 >= m2) {
                merged_system.A(m2, m1) = -2/R;
            } else {
                merged_system.A(m1, m2) = -2/R;
            }
        }

        /*
        // method 2 laplacian with harmonic weights (not working)
        for (auto& interconnection : interconnection_indices_list) {
            uint64_t m1 = interconnection.node1_id + offsets[interconnection.mesh1_id];
            uint64_t m2 = interconnection.node2_id + offsets[interconnection.mesh2_id];
            double R = interconnection.R;

            // find all the adjacent elements to node 1 and node 2
            std::unordered_set<uint64_t> triangles_with_node1;
            std::unordered_set<const CurrentDensityProp*> props1;
            std::unordered_set<uint64_t> triangles_with_node2;
            std::unordered_set<const CurrentDensityProp*> props2;
            // could reserve the size here
            for (uint64_t i = 0; i < meshes[interconnection.mesh1_id].data.numberoftriangles; i++) {
                if (meshes[interconnection.mesh1_id].data.trianglelist[i][0] == interconnection.node1_id ||
                    meshes[interconnection.mesh1_id].data.trianglelist[i][1] == interconnection.node1_id ||
                    meshes[interconnection.mesh1_id].data.trianglelist[i][2] == interconnection.node1_id) {
                    triangles_with_node1.insert(i);
                    props1.insert(meshes[interconnection.mesh1_id].drawing.getRegionPtrFromId(meshes[interconnection.mesh1_id].data.triangleattributelist[i]));
                }
                if (meshes[interconnection.mesh2_id].data.trianglelist[i][0] == interconnection.node2_id ||
                    meshes[interconnection.mesh2_id].data.trianglelist[i][1] == interconnection.node2_id ||
                    meshes[interconnection.mesh2_id].data.trianglelist[i][2] == interconnection.node2_id) {
                    triangles_with_node2.insert(i);
                    props2.insert(meshes[interconnection.mesh2_id].drawing.getRegionPtrFromId(meshes[interconnection.mesh2_id].data.triangleattributelist[i]));
                }
            }

            // nloginfo("found %d triangles with node 1", triangles_with_node1.size());
            for (auto& triangle : triangles_with_node1) {
                // nloginfo("triangle %d", triangle);
            }
            // nloginfo("found %d triangles with node 2", triangles_with_node2.size());
            for (auto& triangle : triangles_with_node2) {
                // nloginfo("triangle %d", triangle);
            }

            // find all the adjacent nodes to node 1 and node 2
            std::unordered_set<uint64_t> nodes_connected_to_node1;
            for (auto& triangle : triangles_with_node1) {
                nodes_connected_to_node1.insert(meshes[interconnection.mesh1_id].data.trianglelist[triangle][0]);
                nodes_connected_to_node1.insert(meshes[interconnection.mesh1_id].data.trianglelist[triangle][1]);
                nodes_connected_to_node1.insert(meshes[interconnection.mesh1_id].data.trianglelist[triangle][2]);
            }

            std::unordered_set<uint64_t> nodes_connected_to_node2;
            for (auto& triangle : triangles_with_node2) {
                nodes_connected_to_node2.insert(meshes[interconnection.mesh2_id].data.trianglelist[triangle][0]);
                nodes_connected_to_node2.insert(meshes[interconnection.mesh2_id].data.trianglelist[triangle][1]);
                nodes_connected_to_node2.insert(meshes[interconnection.mesh2_id].data.trianglelist[triangle][2]);
            }
            // remove node1 from the set
            nodes_connected_to_node1.erase(interconnection.node1_id);
            // remove node2 from the set
            nodes_connected_to_node2.erase(interconnection.node2_id);

            // nloginfo("found %d nodes connected to node 1", nodes_connected_to_node1.size());
            for (auto& node : nodes_connected_to_node1) {
                // nloginfo("node %d", node);
            }
            // nloginfo("found %d nodes connected to node 2", nodes_connected_to_node2.size());
            for (auto& node : nodes_connected_to_node2) {
                // nloginfo("node %d", node);
            }

            // now remove all old weights from the merged system
            std::vector<uint64_t> elems_to_remove;
            for (auto& elem : merged_system.A.elems) {
                uint64_t m = elem.first >> 32;
                uint64_t n = elem.first & 0xFFFFFFFF;

                if (m == interconnection.node1_id + offsets[interconnection.mesh1_id]) {
                    nloginfo("removed weight at %d, %d", m, n);
                    elems_to_remove.push_back(elem.first);
                } else if (m == interconnection.node2_id + offsets[interconnection.mesh2_id]) {
                    nloginfo("removed weight at %d, %d", m, n);
                    elems_to_remove.push_back(elem.first);
                }
            }

            for (auto& elem : elems_to_remove) {
                merged_system.A.elems.erase(elem);
            }

            nloginfo("removed old weights");

            // now add the new weights using the harmonic laplacian approximation
            auto data = meshes[interconnection.mesh1_id].data;
            double sum = 0;
            for (auto& node : nodes_connected_to_node1) {
                // find all the triangles that contain the node
                uint64_t v1 = m1 - offsets[interconnection.mesh1_id];
                uint64_t v2 = node;

                // find the two triangles that contain the current vertex and the adjacent vertex
                uint32_t t1 = 0, t2 = 0;
                bool found_t1 = false, found_t2 = false;

                for (auto myelem : triangles_with_node1) {
                    if (data.trianglelist[myelem][0] == v1 || data.trianglelist[myelem][1] == v1 || data.trianglelist[myelem][2] == v1) {
                        if (data.trianglelist[myelem][0] == v2 || data.trianglelist[myelem][1] == v2 || data.trianglelist[myelem][2] == v2) {
                            if (!found_t1) {
                                t1 = myelem;
                                found_t1 = true;
                            } else {
                                t2 = myelem;
                                found_t2 = true;
                            }
                        }
                    }
                    if (found_t1 && found_t2) break;
                }

                if (!found_t1 || !found_t2) {
                    // we are at the boundary of the mesh, skip this vertex
                    continue;
                }

                // now we have the two triangles, we can compute the cotangent weights
                uint32_t v3_t1, v3_t2;
                if (data.trianglelist[t1][0] == v1 && data.trianglelist[t1][1] == v2) v3_t1 = data.trianglelist[t1][2];
                else if (data.trianglelist[t1][1] == v1 && data.trianglelist[t1][2] == v2) v3_t1 = data.trianglelist[t1][0];
                else if (data.trianglelist[t1][2] == v1 && data.trianglelist[t1][0] == v2) v3_t1 = data.trianglelist[t1][1];
                else if (data.trianglelist[t1][0] == v2 && data.trianglelist[t1][1] == v1) v3_t1 = data.trianglelist[t1][2];
                else if (data.trianglelist[t1][1] == v2 && data.trianglelist[t1][2] == v1) v3_t1 = data.trianglelist[t1][0];
                else if (data.trianglelist[t1][2] == v2 && data.trianglelist[t1][0] == v1) v3_t1 = data.trianglelist[t1][1];
                else nexit("error: vertex not found in element");

                if (data.trianglelist[t2][0] == v1 && data.trianglelist[t2][1] == v2) v3_t2 = data.trianglelist[t2][2];
                else if (data.trianglelist[t2][1] == v1 && data.trianglelist[t2][2] == v2) v3_t2 = data.trianglelist[t2][0];
                else if (data.trianglelist[t2][2] == v1 && data.trianglelist[t2][0] == v2) v3_t2 = data.trianglelist[t2][1];
                else if (data.trianglelist[t2][0] == v2 && data.trianglelist[t2][1] == v1) v3_t2 = data.trianglelist[t2][2];
                else if (data.trianglelist[t2][1] == v2 && data.trianglelist[t2][2] == v1) v3_t2 = data.trianglelist[t2][0];
                else if (data.trianglelist[t2][2] == v2 && data.trianglelist[t2][0] == v1) v3_t2 = data.trianglelist[t2][1];
                else nexit("error: vertex not found in element");

                Vector p1 = data.pointlist[v1];
                Vector p2 = data.pointlist[v2];
                Vector p3_t1 = data.pointlist[v3_t1];
                Vector p3_t2 = data.pointlist[v3_t2];

                double angle1 = geomAngle(p1, p3_t1, p2);
                double angle2 = geomAngle(p1, p3_t2, p2);
                // nloginfo("angle1: %f degrees", angle1 * 180 / M_PI);
                // nloginfo("angle2: %f degrees", angle2 * 180 / M_PI);

                double w = 0.5 * (cos(angle1) / sin(angle1) + cos(angle2) / sin(angle2));
                // nloginfo("w: %.17g", w);

                // we need sigma that we will approximate with the average of the sigmas
                // of all the adjacent triangles
                double sigma = 0;
                for (auto& prop : props1) {
                    sigma += prop->sigma;
                }
                sigma /= props1.size();
                // nloginfo("sigma: %.17g", sigma);

                // our Resistance constraint can be written as:
                // sigma * R * laplacian(V1) - V1 + V2 = 0
                // the weights are our approximation of the laplacian of V in this node

                w *= sigma * R;
                // we will add the -1 to the diagonal and the +1 to the other node later

                sum += w;

                // add the weight to the merged system
                if (v2 >= v1) {
                    nloginfo("added weight %.17g at %d, %d", w, v1 + offsets[interconnection.mesh1_id], v2 + offsets[interconnection.mesh1_id]);
                    merged_system.A(v1 + offsets[interconnection.mesh1_id], v2 + offsets[interconnection.mesh1_id]) = w;
                }
            }

            // now add the -1 to the diagonal
            nloginfo("added weight %.17g at %d, %d", -sum - 1, m1, m1);
            merged_system.A(m1, m1) = -sum - 1;

            // now add the +1 to the other node
            nloginfo("added weight %.17g at %d, %d", 1.0, m1, m2);
            if (m2 >= m1) merged_system.A(m1, m2) = +1;

            // now do the same for node 2
            data = meshes[interconnection.mesh2_id].data;
            sum = 0;
            for (auto& node : nodes_connected_to_node2) {
                // find all the triangles that contain the node
                uint64_t v1 = m2 - offsets[interconnection.mesh2_id];
                uint64_t v2 = node;

                // find the two triangles that contain the current vertex and the adjacent vertex
                uint32_t t1 = 0, t2 = 0;
                bool found_t1 = false, found_t2 = false;
                for (auto myelem : triangles_with_node2) {
                    if (data.trianglelist[myelem][0] == v1 || data.trianglelist[myelem][1] == v1 || data.trianglelist[myelem][2] == v1) {
                        if (data.trianglelist[myelem][0] == v2 || data.trianglelist[myelem][1] == v2 || data.trianglelist[myelem][2] == v2) {
                            if (!found_t1) {
                                t1 = myelem;
                                found_t1 = true;
                            } else {
                                t2 = myelem;
                                found_t2 = true;
                            }
                        }
                    }
                    if (found_t1 && found_t2) break;
                }

                if (!found_t1 || !found_t2) {
                    // we are at the boundary of the mesh, skip this vertex
                    continue;
                }

                // now we have the two triangles, we can compute the cotangent weights
                uint32_t v3_t1, v3_t2;
                if (data.trianglelist[t1][0] == v1 && data.trianglelist[t1][1] == v2) v3_t1 = data.trianglelist[t1][2];
                else if (data.trianglelist[t1][1] == v1 && data.trianglelist[t1][2] == v2) v3_t1 = data.trianglelist[t1][0];
                else if (data.trianglelist[t1][2] == v1 && data.trianglelist[t1][0] == v2) v3_t1 = data.trianglelist[t1][1];
                else if (data.trianglelist[t1][0] == v2 && data.trianglelist[t1][1] == v1) v3_t1 = data.trianglelist[t1][2];
                else if (data.trianglelist[t1][1] == v2 && data.trianglelist[t1][2] == v1) v3_t1 = data.trianglelist[t1][0];
                else if (data.trianglelist[t1][2] == v2 && data.trianglelist[t1][0] == v1) v3_t1 = data.trianglelist[t1][1];
                else nexit("error: vertex not found in element");

                if (data.trianglelist[t2][0] == v1 && data.trianglelist[t2][1] == v2) v3_t2 = data.trianglelist[t2][2];
                else if (data.trianglelist[t2][1] == v1 && data.trianglelist[t2][2] == v2) v3_t2 = data.trianglelist[t2][0];
                else if (data.trianglelist[t2][2] == v1 && data.trianglelist[t2][0] == v2) v3_t2 = data.trianglelist[t2][1];
                else if (data.trianglelist[t2][0] == v2 && data.trianglelist[t2][1] == v1) v3_t2 = data.trianglelist[t2][2];
                else if (data.trianglelist[t2][1] == v2 && data.trianglelist[t2][2] == v1) v3_t2 = data.trianglelist[t2][0];
                else if (data.trianglelist[t2][2] == v2 && data.trianglelist[t2][0] == v1) v3_t2 = data.trianglelist[t2][1];
                else nexit("error: vertex not found in element");

                Vector p1 = data.pointlist[v1];
                Vector p2 = data.pointlist[v2];
                Vector p3_t1 = data.pointlist[v3_t1];
                Vector p3_t2 = data.pointlist[v3_t2];

                double angle1 = geomAngle(p1, p3_t1, p2);
                double angle2 = geomAngle(p1, p3_t2, p2);
                // nloginfo("angle1: %f degrees", angle1 * 180 / M_PI);
                // nloginfo("angle2: %f degrees", angle2 * 180 / M_PI);

                double w = 0.5 * (cos(angle1) / sin(angle1) + cos(angle2) / sin(angle2));

                // nloginfo("w: %.17g", w);

                // we need sigma that we will approximate with the average of the sigmas
                // of all the adjacent triangles
                double sigma = 0;
                for (auto& prop : props2) {
                    sigma += prop->sigma;
                }
                sigma /= props2.size();
                // nloginfo("sigma: %.17g", sigma);

                // our Resistance constraint can be written as:
                // sigma * R * laplacian(V1) - V1 + V2 = 0
                // the weights are our approximation of the laplacian of V in this node

                w *= sigma * R;
                // we will add the -1 to the diagonal and the +1 to the other node later

                sum += w;

                // add the weights to the matrix
                if (v2 >= v1) {
                    nloginfo("added weight %.17g at %d, %d", w, v1 + offsets[interconnection.mesh2_id], v2 + offsets[interconnection.mesh2_id]);
                    merged_system.A(v1 + offsets[interconnection.mesh2_id], v2 + offsets[interconnection.mesh2_id]) = w;
                }
            }

            // now add the -1 to the diagonal
            nloginfo("added weight %.17g at %d, %d", -sum - 1, m2, m2);
            merged_system.A(m2, m2) = -sum - 1;

            // now add the +1 to the other node
            nloginfo("added weight %.17g at %d, %d", 1.0, m2, m1);
            if (m1 >= m2) merged_system.A(m2, m1) = +1;
        }
        */

        // nloginfo("generated merged system with %d unknowns", merged_system.b.val.size());

        return merged_system;
    }

    void MultiLayerCurrentDensitySimulation::setVoltage(CurrentDensitySystem& system, Vector p, double V, uint64_t layer_id) {
        // find the closest node
        int32_t closest_node = -1;
        double closest_distance = INFINITY;
        auto& mesh = meshes[layer_id];

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& _mesh : meshes) {
            offsets.push_back(offset);
            offset += _mesh.data.numberofpoints;
        }

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

        closest_node += offsets[layer_id]; // add the offset

        // set the voltage
        system.addDirichletBoundaryCondition(closest_node, V);
    }

    void MultiLayerCurrentDensitySimulation::solve(CurrentDensitySystem& system) {
        V = CV(system.b.val.size());

        MatCSRSymmetric FemMat(system.A);

        auto start = std::chrono::high_resolution_clock::now();
        preconditionedSSORConjugateGradientSolver(FemMat, system.b, V, 1.5, 1e-20, 100000);
        // preconditionedJacobiConjugateGradientSolver(FemMat, system.b, V, 1e-6, 100000);
        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("solver took %f ms", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    }

    void MultiLayerCurrentDensitySimulation::Vplot(uint32_t width, uint32_t height) {
        std::vector<cv::Mat> images;
        for (uint64_t i = 0; i < meshes.size(); i++) {
            images.push_back(cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255)));
        }

        // get 1st and 99th percentile of V
        std::vector<double> _V = std::vector<double>(V.val.begin(), V.val.end());
        std::sort(_V.begin(), _V.end());
        double max_V = _V[_V.size() * 0.9];
        double min_V = _V[_V.size() * 0.1];

        max_V = 1;
        min_V = -1;

        nloginfo("max_V: %f", max_V);
        nloginfo("min_V: %f", min_V);

        // print middle left and right
        nloginfo("V at middle left: %f", _V[_V.size() * 0.3]);
        nloginfo("V at middle right: %f", _V[_V.size() * 0.7]);

        for (uint64_t i = 0; i < meshes.size(); i++) {
            VplotRend(&images[i], width, height, i, max_V, min_V);
            cv::flip(images[i], images[i], 0); // flip image horizontally
        }

        // stack the images horizontally and show them
        cv::Mat image;
        cv::hconcat(images, image);
        cv::imshow("V", image);
        cv::waitKey(0);
    }

    void MultiLayerCurrentDensitySimulation::VplotRend(cv::Mat* image, double width, double height, uint64_t layer_id, double max_V, double min_V) {
        auto& mesh = meshes[layer_id];
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

        // precompute useful offsets
        std::vector<uint64_t> offsets;
        uint64_t offset = 0;

        // system.b == mesh.data.numberofpoints == A.n == A.m
        for (auto& mesh : meshes) {
            offsets.push_back(offset);
            offset += mesh.data.numberofpoints;
        }

        // get the range of values for the given layer
        std::vector<float> _V = std::vector<float>(V.val.begin() + offsets[layer_id], V.val.begin() + offsets[layer_id] + meshes[layer_id].data.numberofpoints);

        // draw the points
        for (uint32_t i = 0; i < mesh.data.numberofpoints; i++) {
            Vector p = mesh.data.pointlist[i];

            auto points = std::vector<cv::Point>();

            cv::Scalar c = val2jet(_V[i], min_V, max_V);

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
}