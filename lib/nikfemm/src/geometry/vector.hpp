#ifndef NIK_VECTOR_HPP
#define NIK_VECTOR_HPP

namespace nikfemm {
    struct Vector {
        double x;
        double y;

        Vector(double x, double y);
        Vector();

        // == operator with epsilon
        bool operator==(const Vector& v) const;
        // != operator with epsilon
        bool operator!=(const Vector& v) const;
    };
}

#endif