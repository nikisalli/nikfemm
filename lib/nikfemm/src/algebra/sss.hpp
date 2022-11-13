#ifndef NIK_SSS_HPP
#define NIK_SSS_HPP

#include <cstdint>

#include "simple_vector.hpp"
#include "coo.hpp"

namespace nikfemm {
    struct CV;

    struct MatSSS {
        uint32_t* row_ptr;
        uint32_t* col_ind;
        double* val;
        double* diag;

        uint32_t nnz;
        uint32_t m; // number of rows and columns

        MatSSS(const MatCOO& coo);
        MatSSS(const MatSSS& sss);
        ~MatSSS();

        void printSSS();
        void print();
        void write_to_file(const char *filename);
        static void copy(MatSSS& result, const MatSSS& mat);

        double operator()(const uint32_t i, const uint32_t j) const;
        // double& operator()(uint32_t i, uint32_t j);

        void multSSORPreconditioner(CV& result, const CV& cv, double omega);

        void preconditionedJacobiConjugateGradientSolver(CV& b, CV& x0, double maxError, uint32_t maxIterations);
        void preconditionedSSORConjugateGradientSolver(CV& b, CV& x0, double omega, double maxError, uint32_t maxIterations);
    };
}

#endif