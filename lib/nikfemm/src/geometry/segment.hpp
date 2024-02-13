#ifndef NIK_SEGMENT_HPP
#define NIK_SEGMENT_HPP

#include "vector.hpp"
#include "common.hpp"

namespace nikfemm {
    struct Segment {
        Vector p1;
        Vector p2;

        Segment(Vector p1, Vector p2);
        double length();
        static inline bool segmentsIntersect(Segment s1, Segment s2);
        static bool segmentsIntersect(Vector s1p1, Vector s1p2, Vector s2p1, Vector s2p2);
        static bool pointInSegmentBB(Vector p, Segment s);
        static double pointSegmentDistance(Vector p, Segment s);
        static double pointSegmentDistance(Vector p, Vector p1, Vector p2);

        bool operator==(const Segment& s) const;
        bool operator!=(const Segment& s) const;
    };
}

namespace std {
    template <>
    struct hash<nikfemm::Segment> {
        std::size_t operator()(const nikfemm::Segment& s) const {
            return hash<nikfemm::Vector>()(s.p1) ^ hash<nikfemm::Vector>()(s.p2);
        }
    };
}

#endif