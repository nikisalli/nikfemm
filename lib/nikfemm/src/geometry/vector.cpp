#include <math.h>

#include <constants.hpp>

#include "vector.hpp"

namespace nikfemm {

    Vector::Vector(double x, double y) : Point(x, y) {
        this->x = x;
        this->y = y;
    }

    Vector::Vector() : Point() {
    }

    double Vector::magnitude() {
        return sqrt(x * x + y * y);
    }

    Vector Vector::versor() {
        double mag = magnitude();
        return Vector(x / mag, y / mag);
    }
}