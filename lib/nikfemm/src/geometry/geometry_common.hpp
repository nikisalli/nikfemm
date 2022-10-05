#ifndef NIK_GEOMETRY_COMMON_HPP
#define NIK_GEOMETRY_COMMON_HPP

#include "point.hpp"

namespace nikfemm {
    enum Orientation {
        COLLINEAR,
        CLOCKWISE,
        COUNTERCLOCKWISE
    };

    inline double geomDistance(Point p1, Point p2) {
        return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    }

    inline double geomAngle(Point a, Point b, Point c) {
        // safe angle calculation
        double ax = a.x - b.x;
        double ay = a.y - b.y;
        double bx = c.x - b.x;
        double by = c.y - b.y;
        double dot = ax * bx + ay * by;
        double det = ax * by - ay * bx;
        double angle = fabs(atan2(det, dot));
        return angle;
    }

    inline double geomArea(Point a, Point b, Point c) {
        double ab = geomDistance(a, b);
        double bc = geomDistance(b, c);
        double ac = geomDistance(a, c);
        double s = (ab + bc + ac) / 2;
        return sqrt(s * (s - ab) * (s - bc) * (s - ac));
    }

    inline Orientation geomOrientation(Point p1, Point p2, Point p3) {
        double val = (p2.y - p1.y) * (p3.x - p2.x) -
                     (p2.x - p1.x) * (p3.y - p2.y);

        if (abs(val) < EPSILON) {
            return Orientation::COLLINEAR;
        } else if (val > 0) {
            return Orientation::CLOCKWISE;
        } else {
            return Orientation::COUNTERCLOCKWISE;
        }
    }
}

#endif