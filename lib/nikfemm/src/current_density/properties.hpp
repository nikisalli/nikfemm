#ifndef NIK_CURRENT_DENSITY_PROPERTIES_HPP
#define NIK_CURRENT_DENSITY_PROPERTIES_HPP

#include "../geometry/vector.hpp"

namespace nikfemm {
    struct CurrentDensityProp {
        float sigma; // conductivity 

        bool operator==(const CurrentDensityProp& p) const;
        bool operator!=(const CurrentDensityProp& p) const;
    };

    // default material property
    namespace current_density_materials {
        const double copper = (1.0 / 16.78e-9); // 16.78 nΩ•m
    }
}

template <>
struct std::hash<nikfemm::CurrentDensityProp> {
    inline std::size_t operator()(const nikfemm::CurrentDensityProp& p) const {
        return std::hash<float>()(p.sigma);
    }
};

#endif