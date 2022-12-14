#include "utils.hpp"

namespace nikfemm {
    void nexit(std::string message) {
        throw std::runtime_error(message);
    }

    void nassert(bool condition, std::string message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    cv::Scalar val2jet(float v, float vmin, float vmax) {
        cv::Scalar c = {255, 255, 255, 255}; // white
        float dv;

        if (v < vmin)
            v = vmin;
        if (v > vmax)
            v = vmax;
        dv = vmax - vmin;

        dv = vmax - vmin;

        // opencv is BGR
        if (v < (vmin + 0.25 * dv)) {
            c[2] = 0;
            c[1] = (4 * (v - vmin) / dv) * 255;
        } else if (v < (vmin + 0.5 * dv)) {
            c[2] = 0;
            c[0] = (1 + 4 * (vmin + 0.25 * dv - v) / dv) * 255;
        } else if (v < (vmin + 0.75 * dv)) {
            c[2] = (4 * (v - vmin - 0.5 * dv) / dv) * 255;
            c[0] = 0;
        } else {
            c[1] = (1 + 4 * (vmin + 0.75 * dv - v) / dv) * 255;
            c[0] = 0;
        }

        return (c);
    }
}