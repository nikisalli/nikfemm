#ifndef NIK_MAGNETOSTATIC_PROPERTIES_HPP
#define NIK_MAGNETOSTATIC_PROPERTIES_HPP

#include "../geometry/vector.hpp"

namespace nikfemm {
    struct MagnetostaticProp {
        float mu; // permeability
        float J; // current density
        Vector M; // magnetization
        float A; // magnetic vector potential
        Vector B; // magnetic flux density

        // for sorting
        bool operator==(const MagnetostaticProp& other) const;
        bool operator!=(const MagnetostaticProp& other) const;
        bool operator<(const MagnetostaticProp& other) const;
    };

    // default material property
    const MagnetostaticProp vacuum_prop = {MU_0, 0, Vector(0, 0)};
}

#endif