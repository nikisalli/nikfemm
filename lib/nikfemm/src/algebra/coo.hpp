#ifndef NIK_COO_HPP
#define NIK_COO_HPP

#include <set>
#include <cstdint>
#include <unordered_map>
#include <math.h>

#include "assert.h"

#include "SDL2/SDL.h"

#include "coo.hpp"
#include "../utils/utils.hpp"
namespace nikfemm {
    template <typename T>
    struct MatCOO {
        // store only non-zero elements of symmetric matrices in upper triangular part
        std::unordered_map<uint64_t, T> elems;
        
        uint32_t m = 0;  // rows, columns

        MatCOO(uint32_t m);
        ~MatCOO();

        static inline uint64_t get_key(uint32_t i, uint32_t j) {
            // swap i, j if i > j
            if (i > j) {
                std::swap(i, j);
            }
            return (uint64_t) i << 32 | (uint64_t)j;
        }

        double operator()(uint32_t _m, uint32_t _n) const;
    };

    template <typename T>
    MatCOO<T>::MatCOO(uint32_t m) {
        this->m = m;
    }

    template <typename T>
    MatCOO<T>::~MatCOO() {
    }

    template <typename T>
    double MatCOO<T>::operator()(uint32_t _m, uint32_t _n) const {
        // swap indices if _m > _n
        if (_m > _n) {
            std::swap(_m, _n);
        }
        uint64_t key = (uint64_t)_m << 32 | (uint64_t)_n;
        if (elems.find(key) != elems.end()) {
            return elems.at(key);
        } else {
            return 0.0;
        }
    }
}

#endif