#include "polygon.hpp"
#include "segment.hpp"

#include "../utils/utils.hpp"

namespace nikfemm {
    Polygon::Polygon() {
    }

    Polygon::Polygon(const std::vector<Vector>& points) {
        this->points = points;
    }

    Polygon::Polygon(const Vector* points, size_t n) {
        this->points = std::vector<Vector>(points, points + n);
    }

    Polygon::Polygon(const Polygon& p) {
        points = p.points;
    }

    bool Polygon::operator==(const Polygon& p) const {
        if (points.size() != p.points.size()) {
            return false;
        }

        for (size_t i = 0; i < points.size(); i++) {
            if (points[i] != p.points[i]) {
                return false;
            }
        }

        return true;
    }

    bool Polygon::operator!=(const Polygon& p) const {
        return !(*this == p);
    }

    Polygon Polygon::operator=(const Polygon& p) {
        points = p.points;
        return *this;
    }

    bool Polygon::contains(Vector p) const {
        bool inside = false;

        // Iterate through each edge of the polygon
        for (int i = points.size() - 1, j = 0; j < points.size(); i = j, j++) {
            Vector A = points[i];
            Vector B = points[j];

            // Check for corner cases
            if ((p.x == A.x && p.y == A.y) || (p.x == B.x && p.y == B.y)) return true;
            if (A.y == B.y && p.y == A.y && ((p.x >= A.x && p.x <= B.x) || (p.x <= A.x && p.x >= B.x))) return true;

            // Check if p is within the vertical range of the edge
            if ((p.y >= A.y && p.y <= B.y) || (p.y <= A.y && p.y >= B.y)) {
                // Filter out "ray pass vertex" problem
                if (p.y == A.y && B.y >= A.y || p.y == B.y && A.y >= B.y) continue;

                // Calculate the cross product of vectors PA and PB
                double c = (A.x - p.x) * (B.y - p.y) - (B.x - p.x) * (A.y - p.y);
                if (c == 0) return true;

                // If PA is to the left of AB, toggle the inside flag
                if ((A.y < B.y) == (c > 0)) inside = !inside;
            }
        }

        // Return true if p is inside the polygon, false if it is outside
        return inside;
    }

    bool Polygon::contains(Vector p, bool on_boundary_counts, double epsilon) const {
        bool is_contained_simple = contains(p);
        bool is_on_boundary = false;
        for (size_t i = 0; i < points.size(); i++) {
            Vector A = points[i];
            Vector B = points[(i + 1) % points.size()];

            if (Segment::pointSegmentDistance(p, A, B) <= epsilon) {
                is_on_boundary = true;
                break;
            }
        }
        if (on_boundary_counts) {
            return is_contained_simple || is_on_boundary;
        } else {
            return is_contained_simple && !is_on_boundary;
        }
    }

    bool Polygon::contains(Polygon p) const {
        for (Vector point : p.points) {
            if (!contains(point)) {
                return false;
            }
        }

        return true;
    }

    Polygon& Polygon::translate(const Vector v) {
        for (Vector& point : points) {
            point += v;
        }
        return *this;
    }

    Polygon& Polygon::rotate(const double angle, const Vector center) {
        for (Vector& point : points) {
            point = point.rotate(angle, center);
        }
        return *this;
    }

    void Polygon::print() const {
        nloginfo("Polygon: ");
        for (Vector point : points) {
            printf("(%f, %f) ", point.x, point.y);
        }
        printf("\n");
    }
}