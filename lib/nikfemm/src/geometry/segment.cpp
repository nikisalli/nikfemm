#include <math.h>

#include "segment.hpp"
#include "vector.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    Segment::Segment(Vector p1, Vector p2) {
        this->p1 = p1;
        this->p2 = p2;
    }

    float Segment::length() {
        return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    }

    bool Segment::segmentsIntersect(Segment s1, Segment s2) {
        return segmentsIntersect(s1.p1, s1.p2, s2.p1, s2.p2);
    }

    bool Segment::pointInSegmentBB(Vector p, Segment s) {
        return (p.x <= std::max(s.p1.x, s.p2.x) && p.x >= std::min(s.p1.x, s.p2.x) &&
                p.y <= std::max(s.p1.y, s.p2.y) && p.y >= std::min(s.p1.y, s.p2.y));
    }

    double Segment::pointSegmentDistance(Vector p, Segment s) {
        /*
        double x = p.x;
        double y = p.y;
        double x1 = s.p1.x;
        double y1 = s.p1.y;
        double x2 = s.p2.x;
        double y2 = s.p2.y;

        double A = x - x1;
        double B = y - y1;
        double C = x2 - x1;
        double D = y2 - y1;

        double dot = A * C + B * D;
        double len_sq = C * C + D * D;
        double param = -1;
        if (len_sq != 0) //in case of 0 length line
            param = dot / len_sq;

        double xx, yy;

        if (param < 0) {
            xx = x1;
            yy = y1;
        }
        else if (param > 1) {
            xx = x2;
            yy = y2;
        }
        else {
            xx = x1 + param * C;
            yy = y1 + param * D;
        }

        double dx = x - xx;
        double dy = y - yy;
        return sqrt(dx * dx + dy * dy);
        */
        return pointSegmentDistance(p, s.p1, s.p2);
    }

    double Segment::pointSegmentDistance(Vector p, Vector p1, Vector p2) {
        /*
        // Return minimum distance between line segment vw and point p
        const float l2 = Vector::distance_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
        if (l2 == 0.0) return Vector::distance(p, v);   // v == w case
        // Consider the line extending the segment, parameterized as v + t (w - v).
        // We find projection of point p onto the line. 
        // It falls where t = [(p-v) . (w-v)] / |w-v|^2
        // We clamp t from [0,1] to handle points outside the segment vw.
        const float t = std::max(0.f, std::min(1.f, Vector::dot(p - v, w - v) / l2));
        const Vector projection = v + ((w - v) * t);  // Projection falls on the segment
        return Vector::distance(p, projection);
        */
        // return pointSegmentDistance(p, Segment(p1, p2));
        double param = ((p - p1) * (p2 - p1)) / Vector::distanceSquared(p1, p2);
        if (param < 0) {
            return Vector::distance(p, p1);
        }
        else if (param > 1) {
            return Vector::distance(p, p2);
        }
        else {
            return Vector::distance(p, p1 + ((p2 - p1) * param));
        }
    }

    double area(Vector p1, Vector p2, Vector p3) {
        return (p2.x-p1.x)*(p3.y-p1.y)-(p3.x-p1.x)*(p2.y-p1.y);
    }

    bool Segment::segmentsIntersect(Vector s1p1, Vector s1p2, Vector s2p1, Vector s2p2) {
        double epsilon = 1e-100;
        double c_area = area(s1p1, s1p2, s2p1);
        double d_area = area(s1p1, s1p2, s2p2);
        if (std::abs(c_area) < epsilon) {
            if (std::abs(s2p1.x-s1p1.x) < epsilon) {
                if (std::min(s1p1.y,s1p2.y)-epsilon < s2p1.y && s2p1.y < std::max(s1p1.y,s1p2.y)+epsilon) {
                    return true;
                }
            } else {
                if (std::min(s1p1.x,s1p2.x)-epsilon < s2p1.x && s2p1.x < std::max(s1p1.x,s1p2.x)+epsilon) {
                    return true;
                }
            }
            if (std::abs(d_area) > epsilon) {
                return false;
            }
        }
        if (std::abs(d_area) < epsilon) {
            if (std::abs(s2p2.x-s1p1.x) < epsilon) {
                if (std::min(s1p1.y,s1p2.y)-epsilon < s2p2.y && s2p2.y < std::max(s1p1.y,s1p2.y)+epsilon) {
                    return true;
                }
            } else {
                if (std::min(s1p1.x,s1p2.x)-epsilon < s2p2.x && s2p2.x < std::max(s1p1.x,s1p2.x)+epsilon) {
                    return true;
                }
            }
            if (std::abs(c_area) > epsilon) {
                return false;
            }
            if (std::abs(s2p1.x-s1p1.x) < epsilon) {
                return (s1p1.y < s2p1.y) != (s1p1.y < s2p2.y);
            } else {
                return (s1p1.x < s2p1.x) != (s1p1.x < s2p2.x);
            }
        }
        if ((c_area > 0) == (d_area > 0)) {
            return false;
        }
        double a_area = area(s2p1, s2p2, s1p1);
        double b_area = area(s2p1, s2p2, s1p2);
        return (a_area > 0) != (b_area > 0);
    }

    bool Segment::operator==(const Segment& s) const {
        return (p1 == s.p1 && p2 == s.p2) || (p1 == s.p2 && p2 == s.p1);
    }

    bool Segment::operator!=(const Segment& s) const {
        return !(*this == s);
    }
}