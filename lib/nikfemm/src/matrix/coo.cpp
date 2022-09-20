#include <math.h>

#include "coo.hpp"

namespace nikfemm {
    bool ElemCOO::operator<(const ElemCOO &other) const {
        if (m < other.m) {
            return true;
        } else if (m == other.m) {
            return n < other.n;
        } else {
            return false;
        }
    }

    bool ElemCOO::operator==(const ElemCOO &other) const {
        return m == other.m && n == other.n;
    }

    MatCOO::MatCOO() {
    }

    MatCOO::~MatCOO() {
    }

    void MatCOO::add_elem(uint64_t m, uint64_t n, double val) {
        add_elem(ElemCOO{m, n, val});
    }

    void MatCOO::add_elem(ElemCOO elem) {
        if (elem.m + 1 > m ) {
            m = elem.m + 1;
        }
        if (elem.n + 1 > n) {
            n = elem.n + 1;
        }
        if (elem.val != 0.0) {
            elems.push_back(elem);
        }
    }
}