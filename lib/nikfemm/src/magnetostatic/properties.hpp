#ifndef NIK_MAGNETOSTATIC_PROPERTIES_HPP
#define NIK_MAGNETOSTATIC_PROPERTIES_HPP

#include <set>

#include "../geometry/vector.hpp"

namespace nikfemm {
    struct BH_Point {
        float B;
        float H;

        bool operator==(const BH_Point& p) const;
        bool operator!=(const BH_Point& p) const;
    };
}

template <>
struct std::less<nikfemm::BH_Point> {
    inline bool operator()(const nikfemm::BH_Point& p1, const nikfemm::BH_Point& p2) const {
        return p1.B < p2.B;
    }
};

template <>
struct std::equal_to<nikfemm::BH_Point> {
    inline bool operator()(const nikfemm::BH_Point& p1, const nikfemm::BH_Point& p2) const {
        return p1.B == p2.B;
    }
};

namespace nikfemm {
    typedef std::set<BH_Point> BH_Curve;

    typedef double MaxwellStressTensor[2][2];

    struct MagnetostaticProp {
        float J; // current density
        Vector M; // magnetization
        float mu; // permeability
        BH_Curve bh_curve; // BH curve

        float getMu(float B) const;
        bool isLinear() const;

        bool operator==(const MagnetostaticProp& p) const;
        bool operator!=(const MagnetostaticProp& p) const;
    };

    // default material property
    namespace materials {
        const double vacuum = 4 * PI * 1e-7;
        const double air = 4 * PI * 1e-7;
        const double copper = 1.256629e-6;
        const BH_Curve iron = {{1, 2000}, {1.42, 5000}, {1.7, 10000}, {2, 40000}};
        const double iron_linear = 6.3e-3;
        const double neodymium = 1.32e-6;
    }
}

template <>
struct std::hash<nikfemm::MagnetostaticProp> {
    inline std::size_t operator()(const nikfemm::MagnetostaticProp& p) const {
        return std::hash<float>()(p.mu) ^ std::hash<float>()(p.J) ^ std::hash<nikfemm::Vector>()(p.M);
    }
};

template <>
struct std::equal_to<nikfemm::MagnetostaticProp> {
    inline bool operator()(const nikfemm::MagnetostaticProp& p1, const nikfemm::MagnetostaticProp& p2) const {
        return p1.mu == p2.mu && p1.J == p2.J && p1.M == p2.M && p1.bh_curve == p2.bh_curve;
    }
};

#endif