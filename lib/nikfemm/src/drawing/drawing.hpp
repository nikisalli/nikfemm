#ifndef NIK_DRAWING_HPP
#define NIK_DRAWING_HPP

#include <unordered_set>
#include <utility>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <set>
#include <chrono>
#include <array>

#include "../constants.hpp"

#include "../utils/utils.hpp"
#include "../geometry/vector.hpp"
#include "../geometry/segment.hpp"
#include "drawing_segment.hpp"
#include "../geometry/segment.hpp"
#include "../geometry/circle.hpp"
#include "../geometry/polygon.hpp"

namespace nikfemm {
    typedef std::pair<Vector, uint32_t> DrawingRegion;

    template<typename Prop>
    struct Drawing {
        std::vector<Prop> region_map;
        std::vector<DrawingRegion> regions;
        std::vector<Vector> points;
        std::vector<DrawingSegment> segments;
        std::vector<Polygon> polygons;

        double epsilon;  // The minimum distance between two points of the drawing

        Drawing();
        ~Drawing();

        /* drawing */
        public:
            void drawRectangle(Vector p1, Vector p2);
            void drawRectangle(Vector p1, double width, double height);
            void drawCircle(Vector p, double radius, uint32_t n_segments);
            void drawCircle(Circle c, uint32_t n_segments);
            void drawPolygon(Vector* points, uint32_t n_points);
            void drawPolygon(const std::vector<Vector>& points);
            void drawRegion(Vector p, Prop val);
            const Prop* getRegionPtrFromId(uint32_t id) const;
            Prop getRegionFromId(uint32_t id) const;
            uint32_t getRegionId(Prop val);
            void addRefiningPoints();
        private:
            void drawSegment(Vector p1, Vector p2);
            void drawSegment(Segment s);
    };

    // templated member functions must be defined in the header file
    template <typename Prop>
    Drawing<Prop>::Drawing() {

    }

    template <typename Prop>
    Drawing<Prop>::~Drawing() {
        
    }

    template <typename Prop>
    void Drawing<Prop>::drawRectangle(Vector p1, Vector p2) {
        drawSegment(p1, Vector(p2.x, p1.y));
        drawSegment(Vector(p2.x, p1.y), p2);
        drawSegment(p2, Vector(p1.x, p2.y));
        drawSegment(Vector(p1.x, p2.y), p1);
        polygons.push_back(Polygon({p1, Vector(p2.x, p1.y), p2, Vector(p1.x, p2.y)}));
    }

    template <typename Prop>
    void Drawing<Prop>::drawRectangle(Vector p1, double width, double height) {
        drawRectangle(p1, Vector(p1.x + width, p1.y + height));
    }

    template <typename Prop>
    void Drawing<Prop>::drawCircle(Vector p, double radius, uint32_t n_segments) {
        double angle = 0;
        double angle_step = 2 * PI / n_segments;
        Vector p1 = Vector(p.x + radius, p.y);
        Vector p2;
        std::vector<Vector> points;
        for (uint32_t i = 0; i < n_segments; i++) {
            angle += angle_step;
            p2 = Vector(p.x + radius * cos(angle), p.y + radius * sin(angle));
            drawSegment(p1, p2);
            p1 = p2;
            points.push_back(p2);
        }
        polygons.push_back(Polygon(points));
    }

    template <typename Prop>
    void Drawing<Prop>::drawCircle(Circle c, uint32_t n_segments) {
        drawCircle(c.center, c.radius, n_segments);
    }

    template <typename Prop>
    void Drawing<Prop>::drawPolygon(Vector* points, uint32_t n_points) {
        // check if the polygon self-intersects
        for (uint32_t i = 0; i < n_points; i++) {
            for (uint32_t j = i + 2; j < n_points; j++) {
                if (i == 0 && j == n_points - 1) {
                    continue;
                }
                if (Segment::segmentsIntersect(points[i], points[(i + 1) % n_points], points[j], points[(j + 1) % n_points])) {
                    nexit("Error: polygon self-intersects");
                }
            }
        }

        for (uint32_t i = 0; i < n_points - 1; i++) {
            drawSegment(points[i], points[i + 1]);
        }
        drawSegment(points[n_points - 1], points[0]);
        polygons.push_back(Polygon(points, n_points));
    }

    template <typename Prop>
    void Drawing<Prop>::drawPolygon(const std::vector<Vector>& points) {
        drawPolygon((Vector*) points.data(), points.size());
    }

    template <typename Prop>
    uint32_t Drawing<Prop>::getRegionId(Prop val) {
        for (uint32_t i = 0; i < region_map.size(); i++) {
            if (region_map[i] == val) {
                return i;
            }
        }
        region_map.push_back(val);
        return region_map.size() - 1;
    }

    template <typename Prop>
    Prop Drawing<Prop>::getRegionFromId(uint32_t id) const {
        for (auto it = region_map.begin(); it != region_map.end(); it++) {
            if (it->second == id) {
                return it->first;
            }
        }
        nexit("Error: region id not found");
        // unreachable
        exit(1);
    }

    template <typename Prop>
    const Prop* Drawing<Prop>::getRegionPtrFromId(uint32_t id) const {
        if (id >= region_map.size()) {
            nexit("Error: region id not found");
        }
        
        return &region_map.at(id);
        // unreachable
        exit(1);
    }

    template <typename Prop>
    void Drawing<Prop>::drawRegion(Vector p, Prop val) {
        uint32_t region_id = getRegionId(val);
        regions.push_back(DrawingRegion(p, region_id));
    }

    template <typename Prop>
    void Drawing<Prop>::drawSegment(Vector p1, Vector p2) {
        // check if point is already in points
        bool p1_found = false;
        bool p2_found = false;
        uint32_t p1_id = 0;
        uint32_t p2_id = 0;
        for (uint32_t i = 0; i < points.size(); i++) {
            if (points[i] == p1) {
                p1_found = true;
                p1_id = i;
            }
            if (points[i] == p2) {
                p2_found = true;
                p2_id = i;
            }
        }
        if (!p1_found) {
            p1_id = points.size();
            points.push_back(p1);
        }
        if (!p2_found) {
            p2_id = points.size();
            points.push_back(p2);
        }
        
        for (auto s : segments) {
            if (s.p1 == p1_id && s.p2 == p2_id) {
                // printf("Warning: segment already exists\n");
                return;
            }
            if (Segment::segmentsIntersect(p1, p2, points[s.p1], points[s.p2]) && s.p1 != p1_id && s.p1 != p2_id && s.p2 != p1_id && s.p2 != p2_id) {
                nexit("Error: segment intersects another segment");
            }
        }
        segments.push_back(DrawingSegment(p1_id, p2_id));
    }

    template <typename Prop>
    void Drawing<Prop>::drawSegment(Segment s) {
        drawSegment(s.p1, s.p2);
    }

    template <typename Prop>
    void Drawing<Prop>::addRefiningPoints() {
        // find shortest segment first
        epsilon = std::numeric_limits<double>::max();
        for (auto it = segments.begin(); it != segments.end(); it++) {
            double len= Vector::distance(points[it->p1], points[it->p2]);
            if (len < epsilon) {
                epsilon = len;
            }
        }

        uint32_t n_points = points.size();

        // for each corner, add as many points as possible and make sure they are at least 20 degrees apart
        // points should also stay away from the segments
        for (uint32_t i = 0; i < n_points; i++) {
            // find all segments that contain this point
            std::vector<DrawingSegment> segments_containing_point;
            for (auto it = segments.begin(); it != segments.end(); it++) {
                if (it->p1 == i || it->p2 == i) {
                    segments_containing_point.push_back(*it);
                }
            }

            // check for all combinations of segments if there is one that makes an angle of less than 60 degrees
            for (uint32_t j = 0; j < segments_containing_point.size(); j++) {
                for (uint32_t k = j + 1; k < segments_containing_point.size(); k++) {
                    Vector p1 = points[segments_containing_point[j].p1];
                    Vector p2 = points[segments_containing_point[j].p2];
                    Vector p3 = points[segments_containing_point[k].p1];
                    Vector p4 = points[segments_containing_point[k].p2];
                    
                    double angle;
                    if (p1 == p3) {
                        angle = Vector::absoluteAngle(p2, p1, p4);
                    } else if (p1 == p4) {
                        angle = Vector::absoluteAngle(p2, p1, p3);
                    } else if (p2 == p3) {
                        angle = Vector::absoluteAngle(p1, p2, p4);
                    } else if (p2 == p4) {
                        angle = Vector::absoluteAngle(p1, p2, p3);
                    } else {
                        nexit("Error: segment not found");
                    }
                    if (angle < PI * (120.0 / 180.0)) {
                        // printf("angle is %f, refining\n", angle);
                        goto refine;
                    }            
                }
            }

            continue;

            refine:

            // printf("point %u has %lu segments\n", i, segments_containing_point.size());

            for (uint32_t j = 0; j < 18; j++) {
                double angle = j * M_PI / 9;
                Vector p = Vector(points[i].x + (epsilon / 100) * cos(angle), points[i].y + (epsilon / 100) * sin(angle));
                
                // check distance to segments
                bool too_close = false;
                for (uint32_t k = 0; k < segments_containing_point.size(); k++) {
                    if (Segment::pointSegmentDistance(p, Segment(points[segments_containing_point[k].p1], points[segments_containing_point[k].p2])) < (epsilon / 100) / 2) {
                        too_close = true;
                        break;
                    }
                }
                if (!too_close) {
                    points.push_back(p);
                }
            }
        }
    }
}

#endif