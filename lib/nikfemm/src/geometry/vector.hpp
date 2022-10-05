#ifndef NIK_VECTOR_HPP
#define NIK_VECTOR_HPP

#include "point.hpp"

namespace nikfemm {
    struct Vector : Point {
        Vector(double x, double y);
        Vector();

        double magnitude();
        Vector versor();
    };
}

#endif