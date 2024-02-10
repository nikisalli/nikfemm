#include "mesh.hpp"

namespace nikfemm {
    MagnetostaticMesh::MagnetostaticMesh() {
        default_prop = {0, {0, 0}, static_cast<float>(magnetostatic_materials::air), {}};
    }

    System<NonLinearExpression> MagnetostaticMesh::getFemSystem() {
        System<NonLinearExpression> system(data.numberofpoints);
        // since the stiffness matrix is symmetric, this function only computes the upper triangular part

        auto start = std::chrono::high_resolution_clock::now();

        // COMPUTE FEM WEIGHTS
        auto adjelems_ids = std::vector<std::array<uint32_t, 18>>(data.numberofpoints);
        auto adjelems_props = std::vector<std::array<const MagnetostaticProp*, 18>>(data.numberofpoints);
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

                double double_oriented_area = data.getDoubleOrientedArea(v1, v2, v3);
                
                if (double_oriented_area < 0) {
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
                if (v1 >= i) system.A(i, v1) += NonLinearTerm({double_area * (b1 * b1 + c1 * c1) * 0.5, adjelems_ids[i][j]});
                if (v2 >= i) system.A(i, v2) += NonLinearTerm({double_area * (b2 * b1 + c2 * c1) * 0.5, adjelems_ids[i][j]});
                if (v3 >= i) system.A(i, v3) += NonLinearTerm({double_area * (b3 * b1 + c3 * c1) * 0.5, adjelems_ids[i][j]});

                // set the b vector
                system.b[i] += (double_area * adjelems_props[i][j]->J) / 6;
            }
        }

        // line integral
        // in this section we find all the edges where there is a surface current (i.e. the edges where the M vector changes from element to element)
        // to do this first of all we find all the segments in the drawing that are subject to a surface current and their outward normal
        struct SurfaceCurrentSegment {
            Vector p1, p2; // the two points that define the segment
            double Jm; // surface current density on the segment (M x n) (pointing normal to the simulation plane)
        };

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

        struct Magnet {
            Polygon polygon;
            Vector M;
            std::vector<SurfaceCurrentSegment> segments;
        };

        std::vector<Magnet> magnets;
        for (auto region_poly_array : polygons_that_contain_magnet_region) {
            if (region_poly_array.second.size() == 1) {
                // this is a magnet
                for (auto point : region_poly_array.second[0].points) {
                    magnets.push_back({region_poly_array.second[0], drawing.region_map[region_poly_array.first].M, {}});
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
                        magnets.push_back({polygon, drawing.region_map[region_poly_array.first].M, {}});
                    }
                }
            }
        }

        // now we have to find the segments that are subject to a surface current and their outward normal
        for (auto &magnet : magnets) {
            for (uint32_t i = 0; i < magnet.polygon.points.size(); i++) {
                uint32_t j = (i + 1) % magnet.polygon.points.size();
                Vector p1 = magnet.polygon.points[i];
                Vector p2 = magnet.polygon.points[j];
                Vector n = (p2 - p1).normal().normalize();
                // we need the outward normal
                Vector test = Vector::midPoint(p1, p2) + (n * epsilon * 0.1);  // 0.1 just to make sure that we get no false negatives
                if (magnet.polygon.contains(test)) {
                    n = n * -1;
                }
                // draw arrow from the middle of p1 and p2 in the direction of n
                Vector mid = Vector::midPoint(p1, p2);
                double Jm = magnet.M ^ n;
                magnet.segments.push_back({p1, p2, Jm});
            }
        }

        // now we have to find the edges that are subject to a surface current
        // i.e. the edges that lie on a segment that is subject to a surface current
        // edges should be unique

        // we need an unique list of edges
        struct Edge {
            uint32_t v1, v2;
        };

        struct EdgeHash {
            static uint64_t edge_unique_hash(uint32_t v1, uint32_t v2) {
                return (uint64_t)v1 << 32 | (uint64_t)v2;
            }
            std::size_t operator()(const Edge& edge) const {
                return edge_unique_hash(edge.v1, edge.v2) * (edge.v1 > edge.v2) + edge_unique_hash(edge.v2, edge.v1) * (edge.v1 <= edge.v2);
            }
        };

        struct EdgeEqual {
            bool operator()(const Edge& edge1, const Edge& edge2) const {
                return (edge1.v1 == edge2.v1 && edge1.v2 == edge2.v2) || (edge1.v1 == edge2.v2 && edge1.v2 == edge2.v1);
            }
        };

        std::unordered_set<Edge, EdgeHash, EdgeEqual> edges;
        edges.reserve((data.numberofpoints - 2) * 3);

        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            edges.insert({(uint32_t)data.trianglelist[i][0], (uint32_t)data.trianglelist[i][1]});
            edges.insert({(uint32_t)data.trianglelist[i][1], (uint32_t)data.trianglelist[i][2]});
            edges.insert({(uint32_t)data.trianglelist[i][2], (uint32_t)data.trianglelist[i][0]});
        }

        nloginfo("number of edges: %d, predicted number of edges: %d", edges.size(), (data.numberofpoints - 2) * 3);
        

        for (auto edge : edges) {
            Vector p1 = data.pointlist[edge.v1];
            Vector p2 = data.pointlist[edge.v2];
            for (auto magnet : magnets) {
                for (auto segment : magnet.segments) {
                    double dist1 = Segment::pointSegmentDistance(p1, segment.p1, segment.p2);
                    double dist2 = Segment::pointSegmentDistance(p2, segment.p1, segment.p2);
                    if (dist1 < epsilon * 0.1 && dist2 < epsilon * 0.1) {
                        // this edge is subject to a surface current
                        double length = Vector::distance(p1, p2);
                        system.b[edge.v1] += length * segment.Jm;
                        system.b[edge.v2] += length * segment.Jm;
                    }
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        nloginfo("FEM matrix construction took %d ms", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        

        return system;
    }

    void MagnetostaticMesh::addDirichletInfiniteBoundaryConditions(System<NonLinearExpression>& system) {
        // find three furthest points from the center
        uint32_t p1 = 0;
        uint32_t p2 = 0;
        uint32_t p3 = 0;
        double d1 = 0;
        double d2 = 0;
        double d3 = 0;
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            double dist = Vector::distance(data.pointlist[i], Vector(0, 0));
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

        // make the magnetic vector potential zero on the three points
        system.addDirichletBoundaryCondition(p1, 0);
        system.addDirichletBoundaryCondition(p2, 0);
        system.addDirichletBoundaryCondition(p3, 0);
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
        nloginfo("Found %lu magnets", magnets.size());
        for (auto polygon : magnets) {
            for (uint32_t i = 0; i < polygon.points.size(); i++) {
                Vector p1 = polygon.points[i];
                Vector p2 = polygon.points[(i + 1) % polygon.points.size()];
                
                double multiplier = 10;
                uint32_t n_points = Vector::distance(p1, p2) / (drawing.epsilon / multiplier);
                // multiply by 1.1 to make sure that the points are far enough to not trigger conformity checks
                // Vector normal = Vector::normal(p1, p2).versor() * (drawing.epsilon / multiplier) * 1.1;
                Vector normal = (p2 - p1).normal().normalize() * (drawing.epsilon / multiplier) * 1.1;

                for (uint32_t j = 0; j < n_points; j++) {
                    Vector mypoints[2];
                    mypoints[0] = Vector::lerp(p1, p2, (double)j / n_points) + normal;
                    mypoints[1] = Vector::lerp(p1, p2, (double)j / n_points) - normal;

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
                                    Vector p3 = other_polygon.points[l];
                                    Vector p4 = other_polygon.points[(l + 1) % other_polygon.points.size()];
                                    if (Vector::distance(mypoints[k], p3) < drawing.epsilon / multiplier ||
                                        Vector::distance(mypoints[k], p4) < drawing.epsilon / multiplier) {
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
                                if (Vector::distance(mypoints[k], point) < drawing.epsilon / multiplier) {
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
}