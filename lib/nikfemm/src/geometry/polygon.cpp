#include "polygon.hpp"
#include "segment.hpp"

namespace nikfemm {
    Polygon::Polygon() {
    }

    Polygon::Polygon(const std::vector<Point>& points) {
        this->points = points;
    }

    Polygon::Polygon(const Point* points, size_t n) {
        this->points = std::vector<Point>(points, points + n);
    }

    Polygon::Polygon(const Polygon& p) {
        this->points = p.points;
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

    bool Polygon::contains(Point p) const {
        bool inside = false;

        // Iterate through each edge of the polygon
        for (int i = points.size() - 1, j = 0; j < points.size(); i = j, j++) {
            Point A = points[i];
            Point B = points[j];

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

    bool Polygon::contains(Polygon p) const {
        for (Point point : p.points) {
            if (!contains(point)) {
                return false;
            }
        }

        return true;
    }
}