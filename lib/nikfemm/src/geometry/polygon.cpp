#include "polygon.hpp"

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
}