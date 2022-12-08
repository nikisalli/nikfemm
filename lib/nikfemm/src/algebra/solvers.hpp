#ifndef NIK_SOLVERS_HPP
#define NIK_SOLVERS_HPP

#include <cstdint>
#include <math.h>
#include <stdio.h>

#include "csr.hpp"
#include "simple_vector.hpp"

namespace nikfemm {
    struct MatCSRSymmetric;
    struct CV;

    void multSSORPreconditioner(const MatCSRSymmetric& A, CV& result, const CV& cv, double omega);
    MatCSRLowerTri incompleteCholeskyDecomposition(const MatCSRSymmetric& A);

    void conjugateGradientSolver(const MatCSRSymmetric& A, const CV& b, CV& x0, double maxError, uint32_t maxIterations);
    void preconditionedJacobiConjugateGradientSolver(const MatCSRSymmetric& A, const CV& b, CV& x0, double maxError, uint32_t maxIterations);
    void preconditionedSSORConjugateGradientSolver(const MatCSRSymmetric& A, const CV& b, CV& x0, double omega, double maxError, uint32_t maxIterations);
    void preconditionedIncompleteCholeskyConjugateGradientSolver(const MatCSRSymmetric& A, const CV& b, CV& x0, double maxError, uint32_t maxIterations);
}

#endif