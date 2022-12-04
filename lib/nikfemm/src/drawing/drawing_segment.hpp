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

#endif