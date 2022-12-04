#ifndef NIK_GEOMETRY_COMMON_HPP
#define NIK_GEOMETRY_COMMON_HPP

#include "point.hpp"

namespace nikfemm {
    enum Orientation {
        COLLINEAR,
        CLOCKWISE,
        COUNTERCLOCKWISE
    };

    inline float geomAngle(Point a, Point b, Point c) {
        // safe angle calculation
        float ax = a.x - b.x;
        float ay = a.y - b.y;
        float bx = c.x - b.x;
        float by = c.y - b.y;
        float dot = ax * bx + ay * by;
        float det = ax * by - ay * bx;
        float angle = fabs(atan2(det, dot));
        return angle;
    }

    inline float geomArea(Point a, Point b, Point c) {
        float ab = Point::distance(a, b);
        float bc = Point::distance(b, c);
        float ac = Point::distance(a, c);
        float s = (ab + bc + ac) / 2;
        return sqrt(s * (s - ab) * (s - bc) * (s - ac));
    }

    inline Orientation geomOrientation(Point p1, Point p2, Point p3) {
        float val = (p2.y - p1.y) * (p3.x - p2.x) -
                     (p2.x - p1.x) * (p3.y - p2.y);

        if (abs(val) < std::numeric_limits<float>::epsilon()) {
            return Orientation::COLLINEAR;
        } else if (val > 0) {
            return Orientation::CLOCKWISE;
        } else {
            return Orientation::COUNTERCLOCKWISE;
        }
    }
}

#endif