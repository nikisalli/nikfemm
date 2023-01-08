#ifndef NIK_BUILD_COO_HPP
#define NIK_BUILD_COO_HPP

#include <set>
#include <cstdint>
#include <unordered_map>
#include <math.h>
#include <algorithm>

#include "assert.h"

#include "build_coo.hpp"
#include "../utils/utils.hpp"
namespace nikfemm {
    template <typename T>
    struct BuildMatCOO {
        // store only non-zero elements of symmetric matrices in upper triangular part
        std::unordered_map<uint64_t, T> elems;
        
        uint32_t m = 0;  // rows, columns

        BuildMatCOO(uint32_t m);
        ~BuildMatCOO();

        T& operator()(uint32_t _m, uint32_t _n);
    };

    template <typename T>
    BuildMatCOO<T>::BuildMatCOO(uint32_t m) {
        this->m = m;
    }

    template <typename T>
    BuildMatCOO<T>::~BuildMatCOO() {
    }

    template <typename T>
    T& BuildMatCOO<T>::operator()(uint32_t _m, uint32_t _n) {
        if (_m > _n) {
            nexit("BuildMatCOO: m > n not allowed");
        }
        uint64_t key = (uint64_t)_m << 32 | (uint64_t)_n;
        // check if key exists
        if (elems.find(key) == elems.end()) {
            elems.insert(std::make_pair(key, T()));
        }
        return elems[key];
    }
}

#endif