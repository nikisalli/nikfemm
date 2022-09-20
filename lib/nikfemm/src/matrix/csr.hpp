#ifndef NIK_AIJ_HPP
#define NIK_AIJ_HPP

#include <cstdint>
#include <array>

#include "simple_vector.hpp"
#include "coo.hpp"

namespace nikfemm {
    struct CV;

    struct MatCSR {
        uint64_t* IA;
        uint64_t* JA;
        double* A;

        uint64_t nnz;
        uint64_t m, n;  // I, J  // rows,  cols

        MatCSR(MatCOO& coo);
        ~MatCSR();

        void add(uint64_t row, uint64_t col, double val);
        void printCSR();
        void print();

        double& operator()(uint64_t i, uint64_t j);
        double operator()(uint64_t i, uint64_t j) const;
        MatCSR operator+(const MatCSR& mat) const;
        MatCSR operator-(const MatCSR& mat) const;
        MatCSR operator*(double d) const;
        MatCSR operator/(double d) const;
        MatCSR operator-(double d) const;
        MatCSR operator+(double d) const;
        MatCSR operator=(const MatCSR& mat);
        
        CV operator*(const CV& cv) const;

        CV conjugateGradientSolve(CV& b, CV& x0, double maxError, uint64_t maxIterations);
    };
}

#endif