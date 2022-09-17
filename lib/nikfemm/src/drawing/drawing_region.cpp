#include "drawing_region.hpp"

namespace nikfemm {
    DrawingRegion::DrawingRegion(Point p, uint64_t region_id) {
        this->p = p;
        this->region_id = region_id;
    }

    DrawingRegion::~DrawingRegion() {
        
    }

    bool DrawingRegion::operator==(const DrawingRegion& dr) const {
        return this->p == dr.p && this->region_id == dr.region_id;
    }

    bool DrawingRegion::operator!=(const DrawingRegion& dr) const {
        return !(*this == dr);
    }
}