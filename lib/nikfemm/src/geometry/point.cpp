#include <math.h>
#include <vector>

#include <constants.hpp>

#include "point.hpp"
#include "vector.hpp"

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

    Point Point::operator+(const Point& p) const {
        return Point(x + p.x, y + p.y);
    }

    Point Point::operator-(const Point& p) const {
        return Point(x - p.x, y - p.y);
    }

    Point Point::operator+=(const Point& p) {
        x += p.x;
        y += p.y;
        return *this;
    }

    Point Point::operator-=(const Point& p) {
        x -= p.x;
        y -= p.y;
        return *this;
    }

    Point Point::operator+(const Vector& v) const {
        return Point(x + v.x, y + v.y);
    }

    Point Point::operator-(const Vector& v) const {
        return Point(x - v.x, y - v.y);
    }

    Point Point::operator+=(const Vector& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    Point Point::operator-=(const Vector& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    Point Point::operator*(const double d) const {
        return Point(x * d, y * d);
    }

    Point Point::operator/(const double d) const {
        return Point(x / d, y / d);
    }

    Point::operator Vector() const {
        return Vector(x, y);
    }
}