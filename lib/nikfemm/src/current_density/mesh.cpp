#include "mesh.hpp"

namespace nikfemm {
    CurrentDensityMesh::CurrentDensityMesh(double max_triangle_area) {
        // default material property
        default_prop = {static_cast<float>(current_density_materials::copper)};
        this->max_triangle_area = max_triangle_area;
    }

    CurrentDensityMesh::CurrentDensityMesh() {
        // default material property
        default_prop = {static_cast<float>(current_density_materials::copper)};
        max_triangle_area = 1e-0;
    }

    CurrentDensityMesh::~CurrentDensityMesh() {
        
    }

    CurrentDensitySystem CurrentDensityMesh::getFemSystem() {
        CurrentDensitySystem system = {
            BuildMatCOO<double>(data.numberofpoints),
            CV(data.numberofpoints)
        };

        auto start = std::chrono::high_resolution_clock::now();

        // COMPUTE FEM WEIGHTS
        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(data.numberofpoints);
        auto adjelems_props = std::vector<std::array<const CurrentDensityProp*, 18>>(data.numberofpoints);
        auto adjelems_count = std::vector<uint8_t>(data.numberofpoints, 0);

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = data.trianglelist[i][j];
                adjelems_ids[myid][adjelems_count[myid]] = i;
                adjelems_props[myid][adjelems_count[myid]++] = drawing.getRegionPtrFromId(data.triangleattributelist[i]);
            }
        }

        // surface integral
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            for (uint8_t j = 0; j < adjelems_count[i]; j++) {
                uint32_t v1, v2, v3;
                v1 = i;
                Elem myelem = data.trianglelist[adjelems_ids[i][j]];
                if (i == data.trianglelist[adjelems_ids[i][j]][0]) {
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

                double oriented_area = data.getDoubleOrientedArea(v1, v2, v3);

                if (oriented_area < 0) {
                    std::swap(v2, v3);
                }

                double double_area = data.getDoubleOrientedArea(v1, v2, v3);

                // elements are only added if they are in the upper triangle because the matrix is symmetric and this saves half the memory
                double b1 = (data.pointlist[v2].y - data.pointlist[v3].y) / double_area;
                double c1 = (data.pointlist[v3].x - data.pointlist[v2].x) / double_area;
                double b2 = (data.pointlist[v3].y - data.pointlist[v1].y) / double_area;
                double c2 = (data.pointlist[v1].x - data.pointlist[v3].x) / double_area;
                double b3 = (data.pointlist[v1].y - data.pointlist[v2].y) / double_area;
                double c3 = (data.pointlist[v2].x - data.pointlist[v1].x) / double_area;

                if (v1 >= i) system.A(i, v1) += double_area * (b1 * b1 + c1 * c1) * 0.5 * (-adjelems_props[i][j]->sigma);
                if (v2 >= i) system.A(i, v2) += double_area * (b2 * b1 + c2 * c1) * 0.5 * (-adjelems_props[i][j]->sigma);
                if (v3 >= i) system.A(i, v3) += double_area * (b3 * b1 + c3 * c1) * 0.5 * (-adjelems_props[i][j]->sigma);

                // set the b vector
            }
            system.b.val[i] = 0;
        }

        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("FEM matrix construction took %d ms", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        return system;
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

    CurrentDensitySystem CurrentDensityMesh::getFemSystemCotangentWeights() {
        CurrentDensitySystem system = {
            BuildMatCOO<double>(data.numberofpoints),
            CV(data.numberofpoints)
        };

        auto start = std::chrono::high_resolution_clock::now();

        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(data.numberofpoints);
        auto adjelems_props = std::vector<std::array<const CurrentDensityProp*, 18>>(data.numberofpoints);
        auto adjelems_count = std::vector<uint8_t>(data.numberofpoints, 0);

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = data.trianglelist[i][j];
                adjelems_ids[myid][adjelems_count[myid]] = i;
                adjelems_props[myid][adjelems_count[myid]++] = drawing.getRegionPtrFromId(data.triangleattributelist[i]);
            }
        }

        auto adjverts_ids = std::vector<std::array<uint32_t, 18>>(data.numberofpoints);
        auto adjverts_count = std::vector<uint8_t>(data.numberofpoints, 0);

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            for (uint8_t j = 0; j < 3; j++) {
                uint32_t myid = data.trianglelist[i][j];
                adjverts_ids[myid][adjverts_count[myid]++] = data.trianglelist[i][(j + 1) % 3];
                adjverts_ids[myid][adjverts_count[myid]++] = data.trianglelist[i][(j + 2) % 3];
            }
        }

        // COMPUTE FEM WEIGHTS
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            // for each vertex we compute the cotangent weights
            // to do this we iterate over all adjacent vertices and for each we try to find the two adjacent triangles
            // the two adjacent triangles are the ones that both contain the current vertex and the adjacent vertex
            double sum = 0;
            for (uint8_t j = 0; j < adjverts_count[i]; j++) {
                uint32_t v1, v2, v3;
                v1 = i;
                v2 = adjverts_ids[i][j];
                
                // find the two triangles that contain the current vertex and the adjacent vertex
                uint32_t t1 = 0, t2 = 0;
                bool found_t1 = false, found_t2 = false;
                for (uint8_t k = 0; k < adjelems_count[i]; k++) {
                    if (data.trianglelist[adjelems_ids[i][k]][0] == v1 || data.trianglelist[adjelems_ids[i][k]][1] == v1 || data.trianglelist[adjelems_ids[i][k]][2] == v1) {
                        if (data.trianglelist[adjelems_ids[i][k]][0] == v2 || data.trianglelist[adjelems_ids[i][k]][1] == v2 || data.trianglelist[adjelems_ids[i][k]][2] == v2) {
                            if (!found_t1) {
                                t1 = adjelems_ids[i][k];
                                found_t1 = true;
                            } else {
                                t2 = adjelems_ids[i][k];
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

                double w = 0.5 * (cos(angle1) / sin(angle1) + cos(angle2) / sin(angle2));
                sum += w;

                // add the weights to the matrix
                if (v2 >= v1) system.A(v1, v2) += w;
            }

            system.A(i, i) = -sum;
            system.b.val[i] = 0;
        }

        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("FEM matrix construction took %d ms", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

        return system;
    }
}