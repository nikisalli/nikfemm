#ifndef NIK_UTILS_HPP
#define NIK_UTILS_HPP

#include <string>   
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "../geometry/vector.hpp"

namespace nikfemm {
    void nexit(std::string message);
    void nassert(bool condition, std::string message);
    cv::Scalar val2jet(float v, float vmin, float vmax);
}

#endif