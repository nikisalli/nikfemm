#ifndef NIK_DRAWING_OBJECTS_HPP
#define NIK_DRAWING_OBJECTS_HPP

#include <cstdint>
#include <vector>
#include <unordered_set>

#include "../geometry/point.hpp"
#include "../geometry/segment.hpp"

namespace nikfemm {
    struct DrawingRegion {
        Point p;
        uint64_t region_attribute;

        DrawingRegion(Point p, uint64_t region_attribute);
        ~DrawingRegion();

        bool operator==(const DrawingRegion& dr) const;
        bool operator!=(const DrawingRegion& dr) const;
    };
}

namespace std {
    template <>
    struct hash<nikfemm::DrawingRegion> {
        std::size_t operator()(const nikfemm::DrawingRegion& dr) const {
            return hash<nikfemm::Point>()(dr.p) ^ hash<uint64_t>()(dr.region_attribute);
        }
    };
}

#endif