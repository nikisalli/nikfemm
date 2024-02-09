#ifndef NIK_VECTOR_HPP
#define NIK_VECTOR_HPP

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "../constants.hpp"

namespace nikfemm {
    struct Vector;
    class Vector {
        public:
            double x = 0;
            double y = 0;

            Vector(double x, double y);
            Vector();
            ~Vector();

            // operators
            inline bool operator==(const Vector p) const {
                return (abs(x - p.x) < EPSILON) && (abs(y - p.y) < EPSILON);
            }

            inline bool operator!=(const Vector p) const {
                return !(*this == p);
            }

            inline Vector operator+(const Vector p) const {
                return Vector(x + p.x, y + p.y);
            }

            inline Vector operator-(const Vector p) const {
                return Vector(x - p.x, y - p.y);
            }

            inline Vector operator+=(const Vector p) {
                x += p.x;
                y += p.y;
                return *this;
            }

            inline Vector operator-=(const Vector p) {
                x -= p.x;
                y -= p.y;
                return *this;
            }

            inline Vector operator*(const double d) const {
                return Vector(x * d, y * d);
            }

            inline Vector operator/(const double d) const {
                return Vector(x / d, y / d);
            }

            inline Vector operator*=(const double d) {
                x *= d;
                y *= d;
                return *this;
            }

            inline Vector operator/=(const double d) {
                x /= d;
                y /= d;
                return *this;
            }

            inline Vector operator-() const {
                return Vector(-x, -y);
            }

            // dot product
            inline double operator*(const Vector p) const {
                return x * p.x + y * p.y;
            }

            // since the cross product of two vectors is perpendicular to both of them and we are in 2D
            // we can just return the z component of the cross product
            inline double operator^(const Vector p) const {
                return x * p.y - y * p.x;
            }

            // methods
            inline double norm() const {
                return sqrt(x * x + y * y);
            }

            inline double normSquared() const {
                return x * x + y * y;
            }

            inline Vector normalize() {
                double len = norm();
                if (len == 0) {
                    return *this;
                }
                x /= len;
                y /= len;
                return *this;
            }

            inline Vector rotate(double angle) const {
                double s = sin(angle);
                double c = cos(angle);
                return Vector(x * c - y * s, x * s + y * c);
            }

            inline Vector rotate(double angle, Vector center) const {
                return (*this - center).rotate(angle) + center;
            }

            inline Vector normal() const {
                return Vector(-y, x);
            }

            static inline double doubleOrientedArea(const Vector p1, const Vector p2, const Vector p3) {
                return p1.x * p2.y - p1.x * p3.y - p2.x * p1.y + p2.x * p3.y + p3.x * p1.y - p3.x * p2.y;
            }

            static inline double orientedArea(const Vector p1, const Vector p2, const Vector p3) {
                return doubleOrientedArea(p1, p2, p3) / 2;
            }

            static inline double area(const Vector p1, const Vector p2, const Vector p3) {
                return fabs(orientedArea(p1, p2, p3));
            }

            static inline double distance(const Vector p1, const Vector p2) {
                return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
            }

            static inline double distanceSquared(const Vector p1, const Vector p2) {
                return pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
            }

            static inline Vector lerp(const Vector p1, const Vector p2, const double t) {
                return Vector(p1.x + (p2.x - p1.x) * t, p1.y + (p2.y - p1.y) * t);
            }

            static inline Vector midPoint(const Vector p1, const Vector p2) {
                return Vector((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
            }

            static inline float absoluteAngle(Vector a, Vector b, Vector c) {
                return fabs(orientedAngle(a, b, c));
            }

            static inline double orientedAngle(const Vector a, const Vector b, const Vector c) {
                // safe angle calculation
                float ax = a.x - b.x;
                float ay = a.y - b.y;
                float bx = c.x - b.x;
                float by = c.y - b.y;
                float dot = ax * bx + ay * by;
                float det = ax * by - ay * bx;
                float angle = atan2(det, dot);
                return angle;
            }
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Vector> {
        inline std::size_t operator()(const nikfemm::Vector& p) const {
            return hash<double>()(p.x) ^ hash<double>()(p.y);
        }
    };
}

#endif