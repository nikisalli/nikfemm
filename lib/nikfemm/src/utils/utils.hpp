#ifndef NIK_UTILS_HPP
#define NIK_UTILS_HPP

#include <string>   
#include <iostream>

#ifdef NIKFEMM_USE_OPENCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#endif

#include "../geometry/vector.hpp"

namespace nikfemm {
    void nexit(std::string message);
    void nassert(bool condition, std::string message);
    void nloginfo(std::string format, ...);
    void nlogwarn(std::string format, ...);
    void nlogerror(std::string format, ...);

#ifdef NIKFEMM_USE_OPENCV
    cv::Scalar val2jet(float v, float vmin, float vmax);
#endif

    inline double map(double x, double in_min, double in_max, double out_min, double out_max) {
        return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }
    inline double limit(double x, double min, double max) {
        if (x < min) {
            return min;
        } else if (x > max) {
            return max;
        } else {
            return x;
        }
    }
}

#endif