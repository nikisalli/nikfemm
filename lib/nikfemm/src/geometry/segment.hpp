#ifndef NIK_SEGMENT_HPP
#define NIK_SEGMENT_HPP

#include "point.hpp"
#include "geometry_common.hpp"

namespace nikfemm {
    struct Segment {
        Point p1;
        Point p2;

        Segment(Point p1, Point p2);
        float length();
        static inline bool segmentsIntersect(Segment s1, Segment s2);
        static bool segmentsIntersect(Point s1p1, Point s1p2, Point s2p1, Point s2p2);
        static bool pointOnSegment(Point p, Segment s);
        static double pointSegmentDistance(Point p, Segment s);

        bool operator==(const Segment& s) const;
        bool operator!=(const Segment& s) const;
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Segment> {
        std::size_t operator()(const nikfemm::Segment& s) const {
            return hash<nikfemm::Point>()(s.p1) ^ hash<nikfemm::Point>()(s.p2);
        }
    };
}

#endif