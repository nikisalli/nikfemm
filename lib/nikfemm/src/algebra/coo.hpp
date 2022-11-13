#ifndef NIK_COO_HPP
#define NIK_COO_HPP

#include <set>
#include <cstdint>
#include <unordered_map>

namespace nikfemm {
    struct MatCOO {
        // store only non-zero elements of symmetric matrices in upper triangular part
        std::unordered_map<uint64_t, double> elems;
        
        uint32_t m = 0;  // rows, columns

        MatCOO(uint32_t m);
        ~MatCOO();

        double operator()(uint32_t _m, uint32_t _n) const;

        void set_elem(uint32_t _m, uint32_t _n, double val);  // row, column, value
        void add_elem(uint32_t _m, uint32_t _n, double val);  // row, column, value
        double get_elem(uint32_t _m, uint32_t _n);  // row, column, value
        void plot();
        void print();
    };
}

#endif