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

template <>
struct std::hash<nikfemm::Vector> {
    inline std::size_t operator()(const nikfemm::Vector& v) const {
        return std::hash<float>()(v.x) ^ std::hash<float>()(v.y);
    }
};

#endif