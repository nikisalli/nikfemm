#include <math.h>

#include <constants.hpp>

#include "vector.hpp"

namespace nikfemm {

    Vector::Vector(float x, float y) : Point(x, y) {
        this->x = x;
        this->y = y;
    }

    Vector::Vector() : Point() {
    }

    float Vector::magnitude() {
        return sqrt(x * x + y * y);
    }

    Vector Vector::versor() {
        float mag = magnitude();
        return Vector(x / mag, y / mag);
    }
}