#include "drawing_region.hpp"

namespace nikfemm {
    DrawingRegion::DrawingRegion(Point p, uint64_t region_attribute) {
        this->p = p;
        this->region_attribute = region_attribute;
    }

    DrawingRegion::~DrawingRegion() {
        
    }

    bool DrawingRegion::operator==(const DrawingRegion& dr) const {
        return this->p == dr.p && this->region_attribute == dr.region_attribute;
    }

    bool DrawingRegion::operator!=(const DrawingRegion& dr) const {
        return !(*this == dr);
    }
}