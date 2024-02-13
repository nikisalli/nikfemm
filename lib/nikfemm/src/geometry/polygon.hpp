#ifndef NIK_POLYGON_HPP
#define NIK_POLYGON_HPP

#include <vector>

#include "vector.hpp"

namespace nikfemm {
    class Polygon {
    public:
        std::vector<Vector> points;

        Polygon();
        Polygon(const std::vector<Vector>& points);
        Polygon(const Vector* points, size_t n);
        Polygon(const Polygon& p);

        bool operator==(const Polygon& p) const;
        bool operator!=(const Polygon& p) const;
        Polygon operator=(const Polygon& p);

        bool contains(Vector p) const;
        bool contains(Vector p, bool on_boundary_counts, double epsilon) const;
        bool contains(Polygon p) const;

        Polygon& translate(const Vector v);
        Polygon& rotate(const double angle, const Vector center);

        void print() const;

        inline Polygon copy() const {
            return Polygon(points);
        }
    };
}

#endif