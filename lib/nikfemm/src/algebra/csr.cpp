#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <math.h>

#include <constants.hpp>

#include "csr.hpp"

#include "../utils/utils.hpp"

namespace nikfemm {
    MatCSR::MatCSR(MatCOO& coo) {
        m = coo.m;

        std::vector<std::pair<uint64_t, double>> elems;
        elems.reserve(coo.elems.size() - m);

        for (auto const& [key, value] : coo.elems) {
            uint64_t _m = key >> 32;
            uint64_t _n = key & 0xFFFFFFFF;
            if (_m < _n) {  // lower triangle
                elems.push_back(std::make_pair(key, value));
                elems.push_back(std::make_pair((_n << 32) | _m, value));  // make symmetric
            } else if (_m == _n) {
                elems.push_back(std::make_pair(key, value));
            }
        }

        // sort elems by key
        std::sort(elems.begin(), elems.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        nnz = elems.size();

        row_ptr = new uint32_t[m + 1]();
        col_ind = new uint32_t[nnz];
        val = new double[nnz];
        diag = new double[m];

        uint32_t i = 0;
        for (auto const& [key, value] : elems) {
            uint32_t _m = key >> 32;
            uint32_t _n = key & 0xFFFFFFFF;
            if (_m == _n) {
                diag[_m] = value;
            }
            col_ind[i] = _n;
            val[i] = value;
            row_ptr[_m + 1]++;
            i++;
        }

        for (uint32_t i = 0; i < m; i++) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    MatCSR::~MatCSR() {
        delete[] row_ptr;
        delete[] col_ind;
        delete[] val;
    }

    void MatCSR::printCSR() {
        printf("m: %lu, nnz: %lu\n", m, nnz);
        printf("row_ptr: ");
        for (uint32_t i = 0; i < m + 1; i++) {
            printf("%lu ", row_ptr[i]);
        }
        printf("\n");
        printf("col_ind: ");
        for (uint32_t i = 0; i < nnz; i++) {
            printf("%lu ", col_ind[i]);
        }
        printf("\n");
        printf("A: ");
        for (uint32_t i = 0; i < nnz; i++) {
            printf("%.1f ", val[i]);
        }
        printf("\n");
    }

    void MatCSR::print() {
        // iterate over CSR elements
        uint32_t idx = 0;
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    double MatCSR::operator()(uint32_t i, uint32_t j) const {
        assert(i <= m && j <= m);
        if (i == j) {
            return diag[i];
        } else {
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_ind[k] == j) {
                    return val[k];
                }
            }
        }
        return 0;
    }

    void MatCSR::write_to_file(const char *filename) {
        FILE *f = fopen(filename, "w");
        if (f == NULL) {
            nexit("Error opening file!\n");
        }

        // iterate over CSR elements
        uint64_t idx = 0;
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = 0; j < m; j++) {
                fprintf(f, "%.17g ", (*this)(i, j));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    void MatCSR::multSSORPreconditioner(CV& result, const CV& cv, double omega) {
        double temp = omega * (2 - omega);
        for (uint32_t i = 0; i < m; i++) {
            result[i] = cv[i] * temp;
        }

        // invert lower triangle
        for (uint32_t i = 0; i < m; i++) {
            result[i] /= diag[i];
            // printf("i: %lu, diag[i]: %.1f, result[i]: %.1f\n", i, diag[i], result[i]);
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_ind[k] > i) {
                    // printf("- k: %lu, col_ind[k]: %lu, val[k]: %.1f\n", k, col_ind[k], val[k]);
                    result[col_ind[k]] -= val[k] * result[i] * omega;
                }
            }
        }


        for (uint32_t i = 0; i < m; i++) {
            result[i] *= diag[i];
        }

        // invert upper triangle
        for (int64_t i = m - 1; i >= 0; i--) {
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_ind[k] > i) {
                    result[i] -= val[k] * result[col_ind[k]] * omega;
                }
            }
            result[i] /= diag[i];
        }
    }

    void MatCSR::conjugateGradientSolver(CV& b, CV& x, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV p(b.m);
        CV Ap(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV::copy(p, r);
        double rTr = CV::squareSum(r);
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
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

    void MatCSR::preconditionedJacobiConjugateGradientSolver(CV& b, CV& x, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV P(b.m);
        for (uint32_t i = 0; i < m; i++) {
            P[i] = 1 / (*this)(i, i);
        }
        CV z(b.m);
        CV::mult(z, P, r);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        double rTzold;
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
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

    void MatCSR::preconditionedSSORConjugateGradientSolver(CV& b, CV& x, double omega, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV z(b.m);
        MatCSR::multSSORPreconditioner(z, r, omega);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        double rTzold;
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
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
            MatCSR::multSSORPreconditioner(z, r, omega);
            double beta = CV::dot(r, z) / rTzold;
            CV::addScaled(p, z, beta, p);
        }
    }
}