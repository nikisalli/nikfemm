#ifndef NIK_CSR_HPP
#define NIK_CSR_HPP

#include <cstdint>

#include "coo.hpp"

namespace nikfemm {
    struct BaseCSR {
        std::vector<uint32_t> row_ptr;
        std::vector<uint32_t> col_ind;
        std::vector<double> val;
        std::vector<double> diag;

        uint32_t nnz;
        uint32_t m;  // i, j  // rows,  columns

        BaseCSR();
        BaseCSR(MatCOOSymmetric<double>& coo);
        BaseCSR(const BaseCSR& csr);
        ~BaseCSR();

        void print();
        void printCSR();
        void write_to_file(const char *filename);
#ifdef NIKFEMM_USE_OPENCV
        void plot(std::string filename);
#endif

        double operator()(uint32_t i, uint32_t j) const;
    };

    struct MatCSRSymmetric;
    struct MatCSRLowerTri;
    struct MatCSRUpperTri;

    struct MatCSRSymmetric : virtual BaseCSR {
        MatCSRSymmetric() : BaseCSR() {}
        MatCSRSymmetric(MatCOOSymmetric<double>& coo) : BaseCSR(coo) {}
        MatCSRSymmetric(const BaseCSR& csr) : BaseCSR(csr) {}
        ~MatCSRSymmetric() {}
    };

    struct MatCSRLowerTri : virtual BaseCSR {
        MatCSRLowerTri() : BaseCSR() {}
        MatCSRLowerTri(MatCOOSymmetric<double>& coo) : BaseCSR(coo) {}
        MatCSRLowerTri(const BaseCSR& csr) : BaseCSR(csr) {}
        ~MatCSRLowerTri() {}
    };

    struct MatCSRUpperTri : virtual BaseCSR {
        MatCSRUpperTri() : BaseCSR() {}
        MatCSRUpperTri(MatCOOSymmetric<double>& coo) : BaseCSR(coo) {}
        MatCSRUpperTri(const BaseCSR& csr) : BaseCSR(csr) {}
        ~MatCSRUpperTri() {}
    };
}

#endif