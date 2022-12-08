#ifndef NIK_COO_HPP
#define NIK_COO_HPP

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

        static inline uint64_t get_key(uint32_t i, uint32_t j);
        void set_elem(uint32_t i, uint32_t j, T val);

        double operator()(uint32_t _m, uint32_t _n) const;
    };

    template <typename T>
    BuildMatCOO<T>::BuildMatCOO(uint32_t m) {
        this->m = m;
    }

    template <typename T>
    BuildMatCOO<T>::~BuildMatCOO() {
    }

    template <typename T>
    inline uint64_t BuildMatCOO<T>::get_key(uint32_t i, uint32_t j) {
        // swap i, j if i > j
        if (i > j) {
            std::swap(i, j);
        }
        return (uint64_t) i << 32 | (uint64_t)j;
    }

    template <typename T>
    void BuildMatCOO<T>::set_elem(uint32_t i, uint32_t j, T val) {
        uint64_t key = get_key(i, j);
        elems[key] = val;
    }

    template <typename T>
    double BuildMatCOO<T>::operator()(uint32_t _m, uint32_t _n) const {
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