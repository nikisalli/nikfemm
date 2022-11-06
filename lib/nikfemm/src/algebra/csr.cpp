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
        n = coo.n;
        nnz = coo.elems.size();

        IA = new uint64_t[m + 1]();
        JA = new uint64_t[nnz];
        A = new double[nnz];

        uint64_t i = 0;
        for (auto const& [key, val] : coo.elems) {
            // printf("key: %lu, m: %lu, n: %lu, val: %g\n", key, key >> 32, key & 0xFFFFFFFF, val);
            JA[i] = key & 0xFFFFFFFF;
            A[i] = val;
            IA[(key >> 32) + 1]++;
            i++;
        }

        for (uint64_t i = 0; i < m; i++) {
            IA[i + 1] += IA[i];
        }
    }

    MatCSR::~MatCSR() {
        delete[] IA;
        delete[] JA;
        delete[] A;
    }

    void MatCSR::printCSR() {
        printf("m: %lu, n: %lu, nnz: %lu\n", m, n, nnz);
        printf("IA: ");
        for (uint64_t i = 0; i < m + 1; i++) {
            printf("%lu ", IA[i]);
        }
        printf("\n");
        printf("JA: ");
        for (uint64_t i = 0; i < nnz; i++) {
            printf("%lu ", JA[i]);
        }
        printf("\n");
        printf("A: ");
        for (uint64_t i = 0; i < nnz; i++) {
            printf("%.1f ", A[i]);
        }
        printf("\n");
    }

    void MatCSR::print() {
        // iterate over CSR elements
        uint64_t idx = 0;
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = 0; j < n; j++) {
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    double MatCSR::operator()(uint64_t i, uint64_t j) const {
        assert(i <= m && j <= n);
        uint64_t row_start = IA[i];
        uint64_t row_end = IA[i + 1];
        for (uint64_t k = row_start; k < row_end; k++) {
            if (JA[k] == j) {
                return A[k];
            }
        }
        return 0;
    }

    CV MatCSR::getInverseDiagonal() const {
        CV cv(m);
        for (uint64_t i = 0; i < m; i++) {
            cv.set_elem(i, 1.0 / (*this)(i, i));
        }
        return cv;
    }

    void MatCSR::conjugateGradientSolve(CV& b, CV& x, double maxError, uint64_t maxIterations) {
        CV r(b.m);
        CV p(b.m);
        CV Ap(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV::copy(p, r);
        double rTr = CV::squareSum(r);
        for (uint64_t i = 0; i < maxIterations; i++) {
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
#ifdef DEBUG_PRINT
            printf("iteration %lu, error: %f\n", i, sqrt(rTrNew));
#endif
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

    void MatCSR::preconditionedConjugateGradientSolve(CV& b, CV& x, double maxError, uint64_t maxIterations) {
        CV r(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV m = getInverseDiagonal();
        CV z(b.m);
        CV::mult(z, m, r);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        double rTzold;
        for (uint64_t i = 0; i < maxIterations; i++) {
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
#ifdef DEBUG_PRINT
            printf("iteration %lu, error: %f\n", i, sqrt(squareError));
#endif
            if (squareError < maxError * maxError) {
#ifdef DEBUG_PRINT
                printf("converged after %lu iterations\n", i);
                printf("x:\n");
                x.print();
#endif
                break;
            }
            CV::mult(z, m, r);
            double beta = CV::dot(r, z) / rTzold;
            CV::addScaled(p, z, beta, p);
        }
    }

    void MatCSR::write_to_file(const char *filename) {
        FILE *f = fopen(filename, "w");
        if (f == NULL) {
            nexit("Error opening file!\n");
        }

        // iterate over CSR elements
        uint64_t idx = 0;
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = 0; j < n; j++) {
                fprintf(f, "%.17g ", (*this)(i, j));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
}