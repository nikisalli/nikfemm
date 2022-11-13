#ifndef NIK_CSR_HPP
#define NIK_CSR_HPP

#include <cstdint>

#include "simple_vector.hpp"
#include "coo.hpp"

namespace nikfemm {
    struct CV;

    struct MatCSR {
        uint32_t* row_ptr;
        uint32_t* col_ind;
        double* val;
        double* diag;

        uint32_t nnz;
        uint32_t m;  // i, j  // rows,  columns

        MatCSR(MatCOO& coo);
        ~MatCSR();

        void printCSR();
        void print();
        void write_to_file(const char *filename);

        double operator()(uint32_t i, uint32_t j) const;
                
        void multSSORPreconditioner(CV& result, const CV& cv, double omega);

        void conjugateGradientSolver(CV& b, CV& x0, double maxError, uint32_t maxIterations);
        void preconditionedJacobiConjugateGradientSolver(CV& b, CV& x0, double maxError, uint32_t maxIterations);
        void preconditionedSSORConjugateGradientSolver(CV& b, CV& x0, double omega, double maxError, uint32_t maxIterations);
    };
}

#endif