#ifndef NIK_DRAWING_SEGMENT_HPP
#define NIK_DRAWING_SEGMENT_HPP

#include <cstdint>

namespace nikfemm {
    struct DrawingSegment {
        uint32_t p1;
        uint32_t p2;

        DrawingSegment(uint32_t p1, uint32_t p2) {
            this->p1 = p1;
            this->p2 = p2;
        }

        bool operator==(const DrawingSegment& s) const {
            return (p1 == s.p1 && p2 == s.p2) || (p1 == s.p2 && p2 == s.p1);
        }
        bool operator!=(const DrawingSegment& s) const {
            return !(*this == s);
        }
    };
}

namespace std {
    template <>
    struct hash<nikfemm::DrawingSegment> {
        inline std::size_t operator()(const nikfemm::DrawingSegment& s) const {
            return hash<uint32_t>()(s.p1) ^ hash<uint32_t>()(s.p2);
        }
    };

    template <>
    struct equal_to<nikfemm::DrawingSegment> {
        inline bool operator()(const nikfemm::DrawingSegment& s1, const nikfemm::DrawingSegment& s2) const {
            return (s1.p1 == s2.p1 && s1.p2 == s2.p2) || (s1.p1 == s2.p2 && s1.p2 == s2.p1);
        }
    };
}

#endif