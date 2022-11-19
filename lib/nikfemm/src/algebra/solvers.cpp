#include "solvers.hpp"

namespace nikfemm {
    void multSSORPreconditioner(MatCSRSymmetric& A, CV& result, const CV& cv, double omega) {
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

    MatCSRLowerTri incompleteCholeskyDecomposition(MatCSRSymmetric& A) {
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

        /*
        for (uint32_t i = R.m; i-- > 0;) {
            Rinvcoo.set_elem(i, i, 1 / R.diag[i]);
            for (int64_t j = R.m - 1; j >= i + 1; j--) {
                double sum = 0;
                for (uint32_t k = R.row_ptr[i]; k < R.row_ptr[i + 1]; k++) {
                    if (R.col_ind[k] <= j && R.col_ind[k] > i) {
                        printf("R.val (%lu, %lu) R.col_ind[k] (%lu, %lu)\n", i, R.col_ind[k], R.col_ind[k], j);
                        sum += R.val[k] * Rinvcoo(j, R.col_ind[k]);
                    }
                }
                if (sum != 0) {
                    printf("setting Rinvcoo(%lu, %lu) to %.4f\n", i, j, -sum / R(i, i));
                    Rinvcoo.set_elem(i, j, -sum / R.diag[i]);
                }
            }
        }

        for (uint32_t i = 0; i < R.m; i++) {
            Rinvcoo.set_elem(i, i, 1 / R.diag[i]);
            for (uint32_t j = 0; j < i; j++) {
                printf("i: %lu, j: %lu\n", i, j);
                double sum = 0;
                for (uint32_t k = j; k < i; k++) {
                    printf("i: %lu, j: %lu, k: %lu, R(i, k): %.4f, Rinvcoo(j, k): %.4f\n", i, j, k, R(i, k), Rinvcoo(j, k));
                    sum += R(i, k) * Rinvcoo(j, k);
                }
                printf("sum: %.4f\n", sum);
                Rinvcoo.set_elem(i, j, -sum / R.diag[i]);
            }
        }

        for (uint32_t i = 0; i < Rt.m; i++) {
            RinvcooT.set_elem(i, i, 1 / Rt.diag[i]);
            for (uint32_t j = 0; j < i; j++) {
                double sum = 0;
                for (uint32_t k = j; k < i; k++) {
                    printf("i: %lu, j: %lu, k: %lu, R(i, k): %.1f, Rinvcoo(j, k): %.1f\n", i, j, k, Rt(i, k), RinvcooT(j, k));
                    sum += Rt(k, i) * RinvcooT(k, j);
                }
                printf("sum: %.4f\n", sum);
                RinvcooT.set_elem(j, i, -sum / Rt.diag[i]);
            }
        }
        */

        // printf("R.inv:\n");
        // R.print();

        return R;
    }

    void conjugateGradientSolver(MatCSRSymmetric& A, CV& b, CV& x, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV p(b.m);
        CV Ap(b.m);
        CV::mult(r, A, x);
        CV::sub(r, b, r);
        CV::copy(p, r);
        double rTr = CV::squareSum(r);
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, A, p);
            double pAp = CV::dot(p, Ap);
            if (fabs(pAp) < std::numeric_limits<double>::epsilon()) {
                pAp = std::numeric_limits<double>::epsilon();
                printf("warning: pAp is zero. approximating with epsilon");
            }
            double alpha = rTr / pAp;
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            double rTrNew = CV::squareSum(r);
            printf("iteration %lu, error: %f\n", i, sqrt(rTrNew));
            if (rTrNew < maxError * maxError) {
            #ifdef DEBUG_PRINT
                printf("converged after %lu iterations\n", i);
                printf("x:\n");
                x.print();
            #endif
                break;
            }
            double beta = rTrNew / rTr;
            CV::addScaled(p, r, beta, p);
            rTr = rTrNew;
        }
    }

    void preconditionedJacobiConjugateGradientSolver(MatCSRSymmetric& A, CV& b, CV& x, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, A, x);
        CV::sub(r, b, r);
        CV P(b.m);
        for (uint32_t i = 0; i < A.m; i++) {
            P[i] = 1 / A.diag[i];
        }
        CV z(b.m);
        CV::mult(z, P, r);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        double rTzold;
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, A, p);
            double alpha = CV::dot(r, z) / CV::dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<double>::epsilon()) {
                alpha = std::numeric_limits<double>::epsilon();
                printf("warning: alpha is zero. approximating with epsilon\n");
            }
            rTzold = CV::dot(r, z);
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            double squareError = CV::squareSum(r);
            printf("iteration %lu, error: %f\n", i, sqrt(squareError));
            if (squareError < maxError * maxError) {
                #ifdef DEBUG_PRINT
                printf("converged after %lu iterations\n", i);
                printf("x:\n");
                x.print();
                #endif
                break;
            }
            CV::mult(z, P, r);
            double beta = CV::dot(r, z) / rTzold;
            CV::addScaled(p, z, beta, p);
        }
    }

    void preconditionedSSORConjugateGradientSolver(MatCSRSymmetric& A, CV& b, CV& x, double omega, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, A, x);
        CV::sub(r, b, r);
        CV z(b.m);
        multSSORPreconditioner(A, z, r, omega);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        double rTzold;
        double squareError;
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, A, p);
            double alpha = CV::dot(r, z) / CV::dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<double>::epsilon()) {
                alpha = std::numeric_limits<double>::epsilon();
                printf("warning: alpha is zero. approximating with epsilon\n");
            }
            rTzold = CV::dot(r, z);
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            squareError = CV::squareSum(r);
            // printf("iteration %lu, error: %f\n", i, sqrt(squareError));
            if (squareError < maxError * maxError) {
                #ifdef DEBUG_PRINT
                printf("converged after %lu iterations\n", i);
                printf("x:\n");
                x.print();
                #endif
                return;
            }
            multSSORPreconditioner(A, z, r, omega);
            double beta = CV::dot(r, z) / rTzold;
            CV::addScaled(p, z, beta, p);
        }
        // printf("didn't converge, last error: %f\n", sqrt(squareError));
    }

    void preconditionedIncompleteCholeskyConjugateGradientSolver(MatCSRSymmetric& A, CV& b, CV& x, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, A, x);
        CV::sub(r, b, r);
        CV z(b.m);
        CV tmp(b.m);
        MatCSRLowerTri L = incompleteCholeskyDecomposition(A);
        CV::mult(tmp, (MatCSRLowerTri)L, r);
        CV::mult(z, (MatCSRUpperTri)L, tmp);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        double rTzold;
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, A, p);
            double alpha = CV::dot(r, z) / CV::dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<double>::epsilon()) {
                alpha = std::numeric_limits<double>::epsilon();
                printf("warning: alpha is zero. approximating with epsilon\n");
            }
            rTzold = CV::dot(r, z);
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            double squareError = CV::squareSum(r);
            printf("iteration %lu, error: %f\n", i, sqrt(squareError));
            if (squareError < maxError * maxError) {
                printf("converged after %lu iterations\n", i);
                #ifdef DEBUG_PRINT
                printf("x:\n");
                x.print();
                #endif
                break;
            }
            CV::mult(tmp, (MatCSRLowerTri)L, r);
            CV::mult(z, (MatCSRUpperTri)L, tmp);
            double beta = CV::dot(r, z) / rTzold;
            CV::addScaled(p, z, beta, p);
        }
    }
}