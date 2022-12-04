#ifndef NIK_VECTOR_HPP
#define NIK_VECTOR_HPP

#include "point.hpp"

namespace nikfemm {
    struct Vector : Point {
        Vector(float x, float y);
        Vector();

        bool operator==(const Vector& v) const;
        bool operator!=(const Vector& v) const;
        Vector operator+(const Vector& v) const;
        Vector operator-(const Vector& v) const;
        Vector operator*(float f) const;
        Vector operator/(float f) const;
        Vector& operator+=(const Vector& v);
        Vector& operator-=(const Vector& v);
        Vector& operator*=(float f);
        Vector& operator/=(float f);

        inline static float dot(const Vector v1, const Vector v2) {
            return v1.x * v2.x + v1.y * v2.y;
        }
        inline static float cross(const Vector v1, const Vector v2) {
            return v1.x * v2.y - v1.y * v2.x;
        }
        inline float magnitude() const {
            return sqrt(x * x + y * y);
        }
        inline Vector versor() const {
            return Vector(x / magnitude(), y / magnitude());
        }
    };
}
template <>
struct std::hash<nikfemm::Vector> {
    inline std::size_t operator()(const nikfemm::Vector& v) const {
        return std::hash<float>()(v.x) ^ std::hash<float>()(v.y);
    }
};

#endif