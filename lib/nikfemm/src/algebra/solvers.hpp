#ifndef NIK_SOLVERS_HPP
#define NIK_SOLVERS_HPP

#include <cstdint>
#include <math.h>
#include <stdio.h>

#include "csr.hpp"

namespace nikfemm {
    struct MatCSRSymmetric;

    void multSSORPreconditioner(const MatCSRSymmetric& A, std::vector<double>& result, const std::vector<double>& cv, double omega);
    MatCSRLowerTri incompleteCholeskyDecomposition(const MatCSRSymmetric& A);

    void conjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x0, double maxError, uint32_t maxIterations);
    void preconditionedJacobiConjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x0, double maxError, uint32_t maxIterations);
    void preconditionedSSORConjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x0, double omega, double maxError, uint32_t maxIterations);
    void preconditionedIncompleteCholeskyConjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x0, double maxError, uint32_t maxIterations);
}

#endif