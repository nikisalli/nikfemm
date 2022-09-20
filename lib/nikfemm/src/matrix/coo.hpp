#ifndef NIK_COO_HPP
#define NIK_COO_HPP

#include <set>
#include <cstdint>
#include <vector>

namespace nikfemm {
    struct ElemCOO {
        uint64_t m;
        uint64_t n;
        double val;

        bool operator<(const ElemCOO &other) const;
        bool operator==(const ElemCOO &other) const;
    };

    struct MatCOO {
        std::vector<ElemCOO> elems;
        
        uint64_t m = 0;
        uint64_t n = 0;

        MatCOO();
        ~MatCOO();

        void add_elem(uint64_t m, uint64_t n, double val);
        void add_elem(ElemCOO elem);
    };
}

#endif