#include <math.h>
#include <vector>

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

    Vector::~Vector() {
    }

}