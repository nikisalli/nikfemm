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

        IA = new uint32_t[m + 1]();
        JA = new uint32_t[nnz];
        A = new float[nnz];

        uint32_t i = 0;
        for (auto const& [key, val] : coo.elems) {
            // printf("key: %lu, m: %lu, n: %lu, val: %g\n", key, key >> 32, key & 0xFFFFFFFF, val);
            JA[i] = key & 0xFFFFFFFF;
            A[i] = val;
            IA[(key >> 32) + 1]++;
            i++;
        }

        for (uint32_t i = 0; i < m; i++) {
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
        for (uint32_t i = 0; i < m + 1; i++) {
            printf("%lu ", IA[i]);
        }
        printf("\n");
        printf("JA: ");
        for (uint32_t i = 0; i < nnz; i++) {
            printf("%lu ", JA[i]);
        }
        printf("\n");
        printf("A: ");
        for (uint32_t i = 0; i < nnz; i++) {
            printf("%.1f ", A[i]);
        }
        printf("\n");
    }

    void MatCSR::print() {
        // iterate over CSR elements
        uint32_t idx = 0;
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < n; j++) {
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    float MatCSR::operator()(uint32_t i, uint32_t j) const {
        assert(i <= m && j <= n);
        uint32_t row_start = IA[i];
        uint32_t row_end = IA[i + 1];
        for (uint32_t k = row_start; k < row_end; k++) {
            if (JA[k] == j) {
                return A[k];
            }
        }
        return 0;
    }

    CV MatCSR::getInverseDiagonal() const {
        CV cv(m);
        for (uint32_t i = 0; i < m; i++) {
            cv.set_elem(i, 1.0 / (*this)(i, i));
        }
        return cv;
    }

    void MatCSR::conjugateGradientSolve(CV& b, CV& x, float maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV p(b.m);
        CV Ap(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV::copy(p, r);
        float rTr = CV::squareSum(r);
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
            float pAp = CV::dot(p, Ap);
            if (fabs(pAp) < std::numeric_limits<float>::epsilon()) {
                pAp = std::numeric_limits<float>::epsilon();
                printf("warning: pAp is zero. approximating with epsilon");
            }
            float alpha = rTr / pAp;
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            float rTrNew = CV::squareSum(r);
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
            float beta = rTrNew / rTr;
            CV::addScaled(p, r, beta, p);
            rTr = rTrNew;
        }
    }

    void MatCSR::preconditionedConjugateGradientSolve(CV& b, CV& x, float maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV m = getInverseDiagonal();
        CV z(b.m);
        CV::mult(z, m, r);
        CV p(b.m);
        CV::copy(p, z);
        CV Ap(b.m);
        float rTzold;
        for (uint32_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
            float alpha = CV::dot(r, z) / CV::dot(p, Ap);
            if (fabs(alpha) < std::numeric_limits<float>::epsilon()) {
                alpha = std::numeric_limits<float>::epsilon();
                printf("warning: alpha is zero. approximating with epsilon\n");
            }
            rTzold = CV::dot(r, z);
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            float squareError = CV::squareSum(r);
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
            float beta = CV::dot(r, z) / rTzold;
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