#ifndef NIK_CSR_HPP
#define NIK_CSR_HPP

#include <cstdint>

#include "simple_vector.hpp"
#include "coo.hpp"

namespace nikfemm {
    struct CV;

    struct MatCSR {
        uint64_t* IA;
        uint64_t* JA;
        double* A;

        uint64_t nnz;
        uint64_t m, n;  // I, J  // rows,  columns

        MatCSR(MatCOO& coo);
        ~MatCSR();

        void printCSR();
        void print();
        void write_to_file(const char *filename);

        double operator()(uint64_t i, uint64_t j) const;
        
        CV operator*(const CV& cv) const;
        
        CV getInverseDiagonal() const;

        void conjugateGradientSolve(CV& b, CV& x0, double maxError, uint64_t maxIterations);
        void preconditionedConjugateGradientSolve(CV& b, CV& x0, double maxError, uint64_t maxIterations);
    };
}

#endif