#ifndef NIK_BUILD_COO_HPP
#define NIK_BUILD_COO_HPP

#include <set>
#include <cstdint>
#include <unordered_map>
#include <math.h>
#include <algorithm>

#include "assert.h"

#include "coo.hpp"
#include "../utils/utils.hpp"
namespace nikfemm {
    template <typename T>
    struct MatCOOSymmetric {
        // store only non-zero elements of symmetric matrices in upper triangular part
        std::unordered_map<uint64_t, T> elems;
        
        uint32_t m = 0;  // rows, columns

        MatCOOSymmetric(uint32_t m);
        ~MatCOOSymmetric();

        T& operator()(uint32_t _m, uint32_t _n);
    };

    template <typename T>
    MatCOOSymmetric<T>::MatCOOSymmetric(uint32_t m) {
        this->m = m;
    }

    template <typename T>
    MatCOOSymmetric<T>::~MatCOOSymmetric() {
    }

    template <typename T>
    T& MatCOOSymmetric<T>::operator()(uint32_t _m, uint32_t _n) {
        if (_m > _n) {
            nexit("MatCOOSymmetric: m > n not allowed");
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