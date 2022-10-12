#include "properties.hpp"

namespace nikfemm {
    bool MagnetostaticProp::operator==(const MagnetostaticProp& other) const {
        return mu == other.mu && J == other.J && M == other.M && A == other.A && B == other.B;
    }

    // for sorting
    bool MagnetostaticProp::operator!=(const MagnetostaticProp& other) const {
        return !(*this == other);
    }

    // for sorting
    bool MagnetostaticProp::operator<(const MagnetostaticProp& other) const {
        if (mu < other.mu) {
            return true;
        } else if (mu > other.mu) {
            return false;
        } else if (J < other.J) {
            return true;
        } else if (J > other.J) {
            return false;
        } else if (M.x < other.M.x) {
            return true;
        } else if (M.x > other.M.x) {
            return false;
        } else if (M.y < other.M.y) {
            return true;
        } else if (M.y > other.M.y) {
            return false;
        } else if (A < other.A) {
            return true;
        } else if (A > other.A) {
            return false;
        } else if (B.x < other.B.x) {
            return true;
        } else if (B.x > other.B.x) {
            return false;
        } else if (B.y < other.B.y) {
            return true;
        } else if (B.y > other.B.y) {
            return false;
        } else {
            return false;
        }
    }
}