#include "element.hpp"

namespace nikfemm {
    Element::~Element() {
        
    }

    double Element::getArea() {
        return 0;
    }

    bool Element::operator==(const Element& e) const {
        return false;
    }

    bool Element::operator!=(const Element& e) const {
        return true;
    }
}