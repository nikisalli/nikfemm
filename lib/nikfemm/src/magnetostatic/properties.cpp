#include "properties.hpp"
#include "../utils/utils.hpp"

namespace nikfemm {
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
}