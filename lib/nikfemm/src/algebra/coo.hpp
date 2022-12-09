#ifndef NIK_COO_HPP
#define NIK_COO_HPP

#include <cstdint>
#include <vector>

#include "simple_vector.hpp"

namespace nikfemm {
    struct ElemCOO {
        uint32_t row;
        uint32_t col;
        double val;
    };

    struct BaseCOO {
        std::vector<ElemCOO> elems;

        uint32_t m;  // square matrix size

        BaseCOO();
        BaseCOO(BuildMatCOO<double>& coo);
        BaseCOO(const BaseCOO& coo);
        ~BaseCOO();

        void printCOO();
        void print();
        void write_to_file(const char *filename);

        double operator()(uint32_t i, uint32_t j) const;
    };

    struct MatCOOSymmetric;
    struct MatCOOLowerTri;
    struct MatCOOUpperTri;

    struct MatCOOSymmetric : virtual BaseCOO {
        MatCOOSymmetric() : BaseCOO() {}
        MatCOOSymmetric(BuildMatCOO<double>& coo) : BaseCOO(coo) {}
        MatCOOSymmetric(const BaseCOO& coo) : BaseCOO(coo) {}
        ~MatCOOSymmetric() {}
    };

    struct MatCOOLowerTri : virtual BaseCOO {
        MatCOOLowerTri() : BaseCOO() {}
        MatCOOLowerTri(BuildMatCOO<double>& coo) : BaseCOO(coo) {}
        MatCOOLowerTri(const BaseCOO& coo) : BaseCOO(coo) {}
        ~MatCOOLowerTri() {}
    };

    struct MatCOOUpperTri : virtual BaseCOO {
        MatCOOUpperTri() : BaseCOO() {}
        MatCOOUpperTri(BuildMatCOO<double>& coo) : BaseCOO(coo) {}
        MatCOOUpperTri(const BaseCOO& coo) : BaseCOO(coo) {}
        ~MatCOOUpperTri() {}
    };
}

#endif  // NIK_COO_HPP