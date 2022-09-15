#include <math.h>

#include <constants.hpp>

#include "vector.hpp"

namespace nikfemm {
    Vector::Vector(double x, double y) {
        this->x = x;
        this->y = y;
    }

    Vector::Vector() {
        this->x = 0;
        this->y = 0;
    }

    bool Vector::operator==(const Vector& v) const {
        return (abs(x - v.x) < EPSILON) && (abs(y - v.y) < EPSILON);
    }

    bool Vector::operator!=(const Vector& v) const {
        return !(*this == v);
    }
}