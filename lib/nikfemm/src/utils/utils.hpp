#ifndef NIK_UTILS_HPP
#define NIK_UTILS_HPP

#include <string>   
#include <iostream>

namespace nikfemm {
    struct vec2 {
        double x, y;
    };

    struct vec3 {
        double x, y, z;
    };

    void nexit(std::string message);
}

#endif