#ifndef NIK_CSR_HPP
#define NIK_CSR_HPP

#include <cstdint>

#include "simple_vector.hpp"
#include "build_coo.hpp"

namespace nikfemm {
    struct CV;

    struct BaseCSR {
        std::vector<uint32_t> row_ptr;
        std::vector<uint32_t> col_ind;
        std::vector<double> val;
        std::vector<double> diag;

        uint32_t nnz;
        uint32_t m;  // i, j  // rows,  columns

        BaseCSR();
        BaseCSR(BuildMatCOO<double>& coo);
        BaseCSR(const BaseCSR& csr);
        ~BaseCSR();
        void printCSR();
    };

    struct MatCSRSymmetric;
    struct MatCSRLowerTri;
    struct MatCSRUpperTri;

    struct MatCSRSymmetric : virtual BaseCSR {
        MatCSRSymmetric();
        MatCSRSymmetric(BuildMatCOO<double>& coo) : BaseCSR(coo) {}
        MatCSRSymmetric(const BaseCSR& csr) : BaseCSR(csr) {}
        ~MatCSRSymmetric() {}

        void print();
        void write_to_file(const char *filename);

        double operator()(uint32_t i, uint32_t j) const;
    };

    struct MatCSRLowerTri : virtual BaseCSR {
        MatCSRLowerTri();
        MatCSRLowerTri(BuildMatCOO<double>& coo) : BaseCSR(coo) {}
        MatCSRLowerTri(const BaseCSR& csr) : BaseCSR(csr) {}
        ~MatCSRLowerTri() {}

        void print();
        void write_to_file(const char *filename);

        double operator()(uint32_t i, uint32_t j) const;
    };

    struct MatCSRUpperTri : virtual BaseCSR {
        MatCSRUpperTri();
        MatCSRUpperTri(BuildMatCOO<double>& coo) : BaseCSR(coo) {}
        MatCSRUpperTri(const BaseCSR& csr) : BaseCSR(csr) {}
        ~MatCSRUpperTri() {}

        void print();
        void write_to_file(const char *filename);

        double operator()(uint32_t i, uint32_t j) const;
    };
}

#endif