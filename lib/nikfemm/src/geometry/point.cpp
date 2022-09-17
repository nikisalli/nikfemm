#include <math.h>
#include <vector>

#include <constants.hpp>

#include "point.hpp"

namespace nikfemm {
    Point::Point(double x, double y) {
        this->x = x;
        this->y = y;
    }

    Point::Point() {
        this->x = 0;
        this->y = 0;
    }

    Point::~Point() {
    }

    bool Point::operator==(const Point& p) const {
        return (abs(x - p.x) < EPSILON) && (abs(y - p.y) < EPSILON);
    }

    bool Point::operator!=(const Point& p) const {
        return !(*this == p);
    }

    bool Point::operator<(const Point& p) const {
        if (abs(x - p.x) < EPSILON) {
            return y < p.y;
        }
        return x < p.x;
    }

    Point Point::operator+(const Point& p) const {
        return Point(x + p.x, y + p.y);
    }

    Point Point::operator-(const Point& p) const {
        return Point(x - p.x, y - p.y);
    }

    Point Point::operator*(const double& d) const {
        return Point(x * d, y * d);
    }

    Point Point::operator/(const double& d) const {
        return Point(x / d, y / d);
    }

    Orientation Point::orientation(Point p1, Point p2, Point p3) {
        double val = (p2.y - p1.y) * (p3.x - p2.x) -
                     (p2.x - p1.x) * (p3.y - p2.y);

        if (abs(val) < EPSILON) {
            return Orientation::COLLINEAR;
        } else if (val > 0) {
            return Orientation::CLOCKWISE;
        } else {
            return Orientation::COUNTERCLOCKWISE;
        }
    }

    std::vector<Point> Point::getWelzlPoints(std::vector<Point> points) {
        std::vector<Point> welzl_points;
        for (int i = 0; i < points.size(); i++) {
            if (points[i] != *this) {
                welzl_points.push_back(points[i]);
            }
        }
        return welzl_points;
    }
}