#include <math.h>

#include "segment.hpp"
#include "vector.hpp"

namespace nikfemm {
    Segment::Segment(Point p1, Point p2) {
        this->p1 = p1;
        this->p2 = p2;
    }

    float Segment::length() {
        return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    }

    bool Segment::segmentsIntersect(Segment s1, Segment s2) {
        return segmentsIntersect(s1.p1, s1.p2, s2.p1, s2.p2);
    }

    bool Segment::pointInSegmentBB(Point p, Segment s) {
        return (p.x <= std::max(s.p1.x, s.p2.x) && p.x >= std::min(s.p1.x, s.p2.x) &&
                p.y <= std::max(s.p1.y, s.p2.y) && p.y >= std::min(s.p1.y, s.p2.y));
    }

    double Segment::pointSegmentDistance(Point p, Segment s) {
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

    double Segment::pointSegmentDistance(Point p, Point p1, Point p2) {
        /*
        // Return minimum distance between line segment vw and point p
        const float l2 = Point::distance_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
        if (l2 == 0.0) return Point::distance(p, v);   // v == w case
        // Consider the line extending the segment, parameterized as v + t (w - v).
        // We find projection of point p onto the line. 
        // It falls where t = [(p-v) . (w-v)] / |w-v|^2
        // We clamp t from [0,1] to handle points outside the segment vw.
        const float t = std::max(0.f, std::min(1.f, Vector::dot(p - v, w - v) / l2));
        const Vector projection = v + ((w - v) * t);  // Projection falls on the segment
        return Point::distance(p, projection);
        */
        // return pointSegmentDistance(p, Segment(p1, p2));
        double param = Vector::dot(p - p1, p2 - p1) / Point::distance_squared(p1, p2);
        if (param < 0) {
            return Point::distance(p, p1);
        }
        else if (param > 1) {
            return Point::distance(p, p2);
        }
        else {
            return Point::distance(p, p1 + ((p2 - p1) * param));
        }
    }

    bool Segment::segmentsIntersect(Point s1p1, Point s1p2, Point s2p1, Point s2p2) {
        // Find the four orientations needed for general and
        // special cases
        Orientation o1 = geomOrientation(s1p1, s1p2, s2p1);
        Orientation o2 = geomOrientation(s1p1, s1p2, s2p2);
        Orientation o3 = geomOrientation(s2p1, s2p2, s1p1);
        Orientation o4 = geomOrientation(s2p1, s2p2, s1p2);
    
        // General case
        if (o1 != o2 && o3 != o4)
            return true;
    
        // Special Cases
        // s1p1, s1p2 and s2p1 are collinear and s2p1 lies on segment s1p1-s1p2
        if (o1 == COLLINEAR && pointInSegmentBB(s2p1, Segment(s1p1, s1p2))) return true;
    
        // s1p1, s1p2 and s2p2 are collinear and s2p2 lies on segment s1p1-s1p2
        if (o2 == COLLINEAR && pointInSegmentBB(s2p2, Segment(s1p1, s1p2))) return true;
    
        // s2p1, s2p2 and s1p1 are collinear and s1p1 lies on segment s2p1-s2p2
        if (o3 == COLLINEAR && pointInSegmentBB(s1p1, Segment(s2p1, s2p2))) return true;
    
        // s2p1, s2p2 and s1p2 are collinear and s1p2 lies on segment s2p1-s2p2
        if (o4 == COLLINEAR && pointInSegmentBB(s1p2, Segment(s2p1, s2p2))) return true;
    
        return false; // Doesn't fall in any of the above cases
    }

    bool Segment::operator==(const Segment& s) const {
        return (p1 == s.p1 && p2 == s.p2) || (p1 == s.p2 && p2 == s.p1);
    }

    bool Segment::operator!=(const Segment& s) const {
        return !(*this == s);
    }
}