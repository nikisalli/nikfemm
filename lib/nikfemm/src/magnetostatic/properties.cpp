#include "properties.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
    bool BH_Point::operator==(const BH_Point& p) const {
        return this->B == p.B && this->H == p.H;
    }

    bool MagnetostaticProp::isLinear() const {
        return bh_curve.size() == 0;
    }

    float MagnetostaticProp::getMu(float B) const {
        // mu = B / H
        // printf("B = %f\n", B);
        if (bh_curve.empty()) {
            // printf("BH curve is empty, returning mu %.17g\n", mu);
            return mu;
        } else if (bh_curve.size() == 1) {
            nexit("BH curve cannot have only one point");
        }
        // interpolate
        BH_Curve::const_iterator it = bh_curve.lower_bound({B, 0});
        // if B is smaller than the smallest B in the curve
        if (it == bh_curve.begin()) {
            float B1 = it->B;
            float H1 = it->H;
            float B2 = (++it)->B;
            float H2 = it->H;

            // mu = dB/dH
            float m = (B2 - B1) / (H2 - H1);
            // printf("B is smaller than the smallest B in the curve, returning mu %.17g\n", m);
            return m;
        } else if (it == bh_curve.end()) {
            float B1 = (--it)->B;
            float H1 = it->H;
            float B2 = (--it)->B;
            float H2 = it->H;

            // mu = dB/dH
            float m = (B2 - B1) / (H2 - H1);
            // printf("B is larger than the largest B in the curve, returning mu %.17g\n", m);
            return m;
        } else {
            // forward difference
            float B1 = it->B;
            float H1 = it->H;
            float B2 = (++it)->B;
            float H2 = it->H;

            // mu = dB/dH
            float m = (B2 - B1) / (H2 - H1);
            // printf("B is in the middle of the curve, returning mu %.17g\n", m);
            return m;
        }
    }

    bool MagnetostaticProp::operator==(const MagnetostaticProp& p) const {
        return mu == p.mu && J == p.J && M == p.M && A == p.A && B == p.B && bh_curve == p.bh_curve;
    }

    bool MagnetostaticProp::operator!=(const MagnetostaticProp& p) const {
        return !(*this == p);
    }
}