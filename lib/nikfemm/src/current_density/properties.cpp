#include "properties.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    bool CurrentDensityProp::operator==(const CurrentDensityProp& p) const {
        return sigma == p.sigma;
    }

    bool CurrentDensityProp::operator!=(const CurrentDensityProp& p) const {
        return !(*this == p);
    }
}