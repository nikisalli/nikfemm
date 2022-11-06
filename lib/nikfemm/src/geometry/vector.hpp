#ifndef NIK_VECTOR_HPP
#define NIK_VECTOR_HPP

#include "point.hpp"

namespace nikfemm {
    struct Vector : Point {
        Vector(float x, float y);
        Vector();

        float magnitude();
        Vector versor();
    };
}

#endif