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

    bool Vector::operator==(const Vector& v) const {
        return x == v.x && y == v.y;
    }

    bool Vector::operator!=(const Vector& v) const {
        return !(*this == v);
    }

    Vector Vector::operator+(const Vector& v) const {
        return Vector(this->x + v.x, this->y + v.y);
    }

    Vector Vector::operator-(const Vector& v) const {
        return Vector(this->x - v.x, this->y - v.y);
    }

    Vector Vector::operator*(float f) const {
        return Vector(this->x * f, this->y * f);
    }

    Vector Vector::operator/(float f) const {
        return Vector(this->x / f, this->y / f);
    }

    Vector& Vector::operator+=(const Vector& v) {
        this->x += v.x;
        this->y += v.y;
        return *this;
    }

    Vector& Vector::operator-=(const Vector& v) {
        this->x -= v.x;
        this->y -= v.y;
        return *this;
    }

    Vector& Vector::operator*=(float f) {
        this->x *= f;
        this->y *= f;
        return *this;
    }

    Vector& Vector::operator/=(float f) {
        this->x /= f;
        this->y /= f;
        return *this;
    }

    float Vector::magnitude() {
        return sqrt(x * x + y * y);
    }

    Vector Vector::versor() {
        float mag = magnitude();
        return Vector(x / mag, y / mag);
    }
}