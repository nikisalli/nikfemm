#ifndef NIK_SSS_HPP
#define NIK_SSS_HPP

#include <cstdint>

#include "simple_vector.hpp"
#include "coo.hpp"

namespace nikfemm {
    struct CV;

    struct MatSSS {
        uint64_t* row_ptr;
        uint64_t* col_ind;
        double* val;
        double* diag;

        uint64_t nnz;
        uint64_t m, n;  // I, J  // rows,  columns

        MatSSS(MatCOO& coo);
        ~MatSSS();

        void printSSS();
        void print();
        void write_to_file(const char *filename);

        double operator()(uint64_t i, uint64_t j) const;
        
        CV getInverseDiagonal() const;

        void conjugateGradientSolve(CV& b, CV& x0, double maxError, uint64_t maxIterations);
        void preconditionedConjugateGradientSolve(CV& b, CV& x0, double maxError, uint64_t maxIterations);
    };
}

#endif