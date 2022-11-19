#include "properties.hpp"

namespace nikfemm {
    float MagnetostaticProp::getMu(float B) const {
        // mu = B / H
        // printf("B = %f\n", B);
        if (bh_curve.empty()) {
            // printf("BH curve is empty, returning mu %.17g\n", mu);
            return mu;
        }
        // interpolate
        BH_Curve::const_iterator it = bh_curve.lower_bound({B, 0});
        // if B is smaller than the smallest B in the curve
        if (it == bh_curve.begin()) {
            // printf("B is smaller than the smallest B in the curve, using {B: %f, H: %f} mu = %f\n", it->B, it->H, it->B / it->H);
            return it->B / it->H;
        } else if (it == bh_curve.end()) {
            // if B is larger than the largest B in the curve
            --it;
            printf("B is larger than the largest B in the curve, using {B: %f, H: %f} mu = %f\n", it->B, it->H, it->B / it->H);
            return it->B / it->H;
        } else {
            BH_Curve::const_iterator it2 = it;
            --it2;
            float B1 = it2->B;
            float B2 = it->B;
            float H1 = it2->H;
            float H2 = it->H;
            float mu1 = B1 / H1;
            float mu2 = B2 / H2;
            float scale = (B - B1) / (B2 - B1);
            float newmu = mu1 + (mu2 - mu1) * scale;

            printf("B is between {B: %f, H: %f} and {B: %f, H: %f} mu = %f\n", B1, H1, B2, H2, newmu);
            return newmu;
        }
    }
}