#include "solvers.hpp"
#include "math.hpp"

namespace nikfemm {
    void multSSORPreconditioner(const MatCSRSymmetric& A, std::vector<double>& result, const std::vector<double>& cv, double omega) {
        double temp = omega * (2 - omega);
        for (uint32_t i = 0; i < A.m; i++) {
            result[i] = cv[i] * temp;
        }

        // invert lower triangle
        for (uint32_t i = 0; i < A.m; i++) {
            result[i] /= A.diag[i];
            // printf("i: %lu, A.diag[i]: %.1f, result[i]: %.1f\n", i, A.diag[i], result[i]);
            for (uint32_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; k++) {
                // printf("- k: %lu, A.col_ind[k]: %lu, A.val[k]: %.1f\n", k, A.col_ind[k], A.val[k]);
                result[A.col_ind[k]] -= A.val[k] * result[i] * omega;
            }
        }

        for (uint32_t i = 0; i < A.m; i++) {
            result[i] *= A.diag[i];
        }

        // invert upper triangle
        for (int64_t i = A.m - 1; i >= 0; i--) {
            for (uint32_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; k++) {
                result[i] -= A.val[k] * result[A.col_ind[k]] * omega;
            }
            result[i] /= A.diag[i];
        }
    }

    MatCSRLowerTri incompleteCholeskyDecomposition(const MatCSRSymmetric& A) {
        MatCSRUpperTri R(A);

        // printf("R:\n");
        // R.print();
        // do incomplete Cholesky factorization using cholesky-crout algorithm
        for (uint32_t i = 0; i < R.m; i++) {
            R.diag[i] = sqrt(R.diag[i]);
            for (uint32_t k = R.row_ptr[i]; k < R.row_ptr[i + 1]; k++) {
                R.val[k] /= R.diag[i];
            }
            for (uint32_t j = i + 1; j < R.m; j++) {
                R.diag[j] -= R(i, j) * R(i, j);
                for (uint32_t k = R.row_ptr[j]; k < R.row_ptr[j + 1]; k++) {
                    if (R.col_ind[k] >= j) {
                        R.val[k] -= R(i, R.col_ind[k]) * R(i, j);
                    }
                }
            }
        }

        // printf("R_chol:\n");
        // R.print();

        // invert R (lower triangular) using forward substitution in place
        for (uint32_t i = R.m; i-- > 0;) {
            R.diag[i] = 1 / R.diag[i];
            for (int64_t j = R.m - 1; j >= i + 1; j--) {
                double sum = 0;
                for (uint32_t k = R.row_ptr[i]; k < R.row_ptr[i + 1]; k++) {
                    if (R.col_ind[k] <= j) {
                        sum += R.val[k] * R(R.col_ind[k], j);
                    }
                }
                if (sum != 0) {
                    for (uint32_t k = R.row_ptr[i]; k < R.row_ptr[i + 1]; k++) {
                        if (R.col_ind[k] == j) {
                            R.val[k] = -sum * R.diag[i];
                        }
                    }
                }
            }
        }

        // printf("R.inv:\n");
        // R.print();

        return R;
    }

    void conjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x, double maxError, uint32_t maxIterations) {
        std::vector<double> r(b.size());
        std::vector<double> p(b.size());
        std::vector<double> Ap(b.size());
        mult(r, A, x);
        sub(r, b, r);
        copy(p, r);
        double rTr = squareSum(r);
        double squareError = squareSum(r);
        for (uint32_t i = 0; i < maxIterations; i++) {
            mult(Ap, A, p);
            double pAp = dot(p, Ap);
            if (fabs(pAp) < std::numeric_limits<double>::epsilon()) {
                pAp = std::numeric_limits<double>::epsilon();
                nlogwarn("warning: pAp is zero. approximating with epsilon");
            }
            double alpha = rTr / pAp;
            addScaled(x, x, alpha, p);
            addScaled(r, r, -alpha, Ap);
            double rTrNew = squareSum(r);
            if (i % 100 == 0) {
                nloginfo("iteration %u, error: %.17g", i, sqrt(squareError));
            }
            if (rTrNew < maxError * maxError) {
            nloginfo("converged after %lu iterations", i);
                // nloginfo("x:");
                // x.print();
                break;
            }
            double beta = rTrNew / rTr;
            addScaled(p, r, beta, p);
            rTr = rTrNew;
        }
    }

    void preconditionedJacobiConjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x, double maxError, uint32_t maxIterations) {
        std::vector<double> r(b.size());
        mult(r, A, x);
        sub(r, b, r);
        std::vector<double> P(b.size());
        for (uint32_t i = 0; i < A.m; i++) {
            P[i] = 1 / A.diag[i];
        }
        std::vector<double> z(b.size());
        mult(z, P, r);
        std::vector<double> p(b.size());
        copy(p, z);
        std::vector<double> Ap(b.size());
        double rTzold;
        double squareError = squareSum(r);
        for (uint32_t i = 0; i < maxIterations; i++) {
            mult(Ap, A, p);
            double alpha = dot(r, z) / dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<double>::epsilon()) {
                alpha = std::numeric_limits<double>::epsilon();
                nlogwarn("warning: alpha is zero. approximating with epsilon");
            }
            rTzold = dot(r, z);
            addScaled(x, x, alpha, p);
            addScaled(r, r, -alpha, Ap);
            squareError = squareSum(r);
            if (i % 100 == 0) {
                nloginfo("iteration %u, error: %.17g", i, sqrt(squareError));
            }
            if (squareError < maxError * maxError) {
                nloginfo("converged after %lu iterations", i);
                // nloginfo("x:");
                // x.print();
                break;
            }
            mult(z, P, r);
            double beta = dot(r, z) / rTzold;
            addScaled(p, z, beta, p);
        }
    }

    void preconditionedSSORConjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x, double omega, double maxError, uint32_t maxIterations) {
        std::vector<double> r(b.size());
        mult(r, A, x);
        sub(r, b, r);
        std::vector<double> z(b.size());
        multSSORPreconditioner(A, z, r, omega);
        std::vector<double> p(b.size());
        copy(p, z);
        std::vector<double> Ap(b.size());
        double rTzold, rTz;
        double squareError = squareSum(r);
        rTz = dot(r, z);
        if (squareError < maxError * maxError) {
            nloginfo("converged after 0 iterations");
            // nloginfo("x:");
            // x.print();
            return;
        }
        for (uint32_t i = 0; i < maxIterations; i++) {
            mult(Ap, A, p);
            double alpha = rTz / dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<double>::epsilon()) {
                alpha = std::numeric_limits<double>::epsilon();
                nlogwarn("warning: alpha is zero. approximating with epsilon");
            }
            addScaled(x, x, alpha, p);
            addScaled(r, r, -alpha, Ap);
            squareError = squareSum(r);
            if (i % 100 == 0) {
                nloginfo("iteration %u, error: %.17g", i, sqrt(squareError));
            }
            if (squareError < maxError * maxError) {
                nloginfo("converged after %lu iterations", i);
                // nloginfo("x:");
                // x.print();
                return;
            }
            multSSORPreconditioner(A, z, r, omega);
            rTzold = rTz;
            rTz = dot(r, z);
            double beta = rTz / rTzold;
            addScaled(p, z, beta, p);
        }
        nloginfo("didn't converge, last error: %.17g", sqrt(squareError));
    }

    void preconditionedIncompleteCholeskyConjugateGradientSolver(const MatCSRSymmetric& A, const std::vector<double>& b, std::vector<double>& x, double maxError, uint32_t maxIterations) {
        std::vector<double> r(b.size());
        mult(r, A, x);
        sub(r, b, r);
        std::vector<double> z(b.size());
        std::vector<double> tmp(b.size());
        nloginfo("computing incomplete cholesky preconditioner");
        MatCSRLowerTri L = incompleteCholeskyDecomposition(A);
        nloginfo("done computing incomplete cholesky preconditioner");
        mult(tmp, (MatCSRLowerTri)L, r);
        mult(z, (MatCSRUpperTri)L, tmp);
        std::vector<double> p(b.size());
        copy(p, z);
        std::vector<double> Ap(b.size());
        double rTzold;
        double squareError = squareSum(r);
        if (squareError < maxError * maxError) {
            nloginfo("converged after 0 iterations");
            // nloginfo("x:");
            // x.print();
            return;
        }
        for (uint32_t i = 0; i < maxIterations; i++) {
            mult(Ap, A, p);
            double alpha = dot(r, z) / dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<double>::epsilon()) {
                alpha = std::numeric_limits<double>::epsilon();
                nlogwarn("warning: alpha is zero. approximating with epsilon");
            }
            rTzold = dot(r, z);
            addScaled(x, x, alpha, p);
            addScaled(r, r, -alpha, Ap);
            squareError = squareSum(r);
            if (i % 100 == 0) {
                nloginfo("iteration %u, error: %.17g", i, sqrt(squareError));
            }
            if (squareError < maxError * maxError) {
                nloginfo("converged after %u iterations", i);
                // nloginfo("x:");
                // x.print();
                break;
            }
            mult(tmp, (MatCSRLowerTri)L, r);
            mult(z, (MatCSRUpperTri)L, tmp);
            double beta = dot(r, z) / rTzold;
            addScaled(p, z, beta, p);
        }
    }
}