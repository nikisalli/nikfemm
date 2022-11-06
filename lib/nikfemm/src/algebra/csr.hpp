#ifndef NIK_CSR_HPP
#define NIK_CSR_HPP

#include <cstdint>

#include "simple_vector.hpp"
#include "coo.hpp"

namespace nikfemm {
    struct CV;

    struct MatCSR {
        uint32_t* IA;
        uint32_t* JA;
        double* A;

        uint32_t nnz;
        uint32_t m, n;  // I, J  // rows,  columns

        MatCSR(MatCOO& coo);
        ~MatCSR();

        void printCSR();
        void print();
        void write_to_file(const char *filename);

        double operator()(uint32_t i, uint32_t j) const;
        
        CV operator*(const CV& cv) const;
        
        CV getInverseDiagonal() const;

        void conjugateGradientSolve(CV& b, CV& x0, double maxError, uint32_t maxIterations);
        void preconditionedConjugateGradientSolve(CV& b, CV& x0, double maxError, uint32_t maxIterations);
    };
}

#endif