#ifndef NIK_GEOMETRY_COMMON_HPP
#define NIK_GEOMETRY_COMMON_HPP

#include "vector.hpp"

namespace nikfemm {
    enum Orientation {
        COLLINEAR,
        CLOCKWISE,
        COUNTERCLOCKWISE
    };

    inline Orientation geomOrientation(Vector p1, Vector p2, Vector p3) {
        double val = (p2.y - p1.y) * (p3.x - p2.x) -
                     (p2.x - p1.x) * (p3.y - p2.y);

        if (abs(val) < std::numeric_limits<double>::epsilon()) {
            return Orientation::COLLINEAR;
        } else if (val > 0) {
            return Orientation::CLOCKWISE;
        } else {
            return Orientation::COUNTERCLOCKWISE;
        }
    }

    inline bool pointInTriangle(Vector p, Vector p1, Vector p2, Vector p3) {
        Orientation o1 = geomOrientation(p, p1, p2);
        Orientation o2 = geomOrientation(p, p2, p3);
        Orientation o3 = geomOrientation(p, p3, p1);

        return ((o1 == o2) && (o2 == o3));
    }
}

#endif