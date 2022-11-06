#ifndef NIK_COO_HPP
#define NIK_COO_HPP

#include <set>
#include <cstdint>
#include <map>

namespace nikfemm {
    struct MatCOO {
        std::map<uint64_t, double> elems;
        
        uint32_t m = 0;  // rows
        uint32_t n = 0;  // columns

        MatCOO(uint32_t m, uint32_t n);
        ~MatCOO();

        void set_elem(uint32_t _m, uint32_t _n, double val);  // row, column, value
        void add_elem(uint32_t _m, uint32_t _n, double val);  // row, column, value
        double get_elem(uint32_t _m, uint32_t _n);  // row, column, value
        void plot();
    };
}

#endif