#include "utils.hpp"

namespace nikfemm {
    void nexit(std::string message) {
        throw std::runtime_error(message);
    }
}