#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>
#include <assert.h>

#include "csr.hpp"

namespace nikfemm {
    MatCSR::MatCSR(MatCOO& coo) {
        // sort elems
        std::sort(coo.elems.begin(), coo.elems.end());
        m = coo.m;
        n = coo.n;
        nnz = coo.elems.size();
        IA = (uint64_t*) calloc(m + 1, sizeof(uint64_t));
        JA = (uint64_t*) calloc(nnz, sizeof(uint64_t));
        A = (double*) calloc(nnz, sizeof(double));

        IA[0] = 0;
        for (uint64_t i = 0; i < nnz; i++) {
            JA[i] = coo.elems[i].n;
            A[i] = coo.elems[i].val;
            IA[coo.elems[i].m + 1]++;
        }
        for (uint64_t i = 0; i < m; i++) {
            IA[i + 1] += IA[i];
        }
    }

    MatCSR::~MatCSR() {
        free(IA);
        free(JA);
        free(A);
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
                if (idx < nnz && JA[idx] == j) {
                    printf("%.1f ", A[idx]);
                    idx++;
                } else {
                    printf("0 ");
                }
            }
            printf("\n");
        }   
    }

    double& MatCSR::operator()(uint64_t i, uint64_t j) {
        // iterate over CSR elements
        uint64_t idx = IA[i];
        while (idx < IA[i + 1]) {
            if (JA[idx] == j) {
                return A[idx];
            }
            idx++;
        }
        throw std::out_of_range("index out of range");
    }

    double MatCSR::operator()(uint64_t i, uint64_t j) const {
        // iterate over CSR elements
        uint64_t idx = IA[i];
        while (idx < IA[i + 1]) {
            if (JA[idx] == j) {
                return A[idx];
            }
            idx++;
        }
        throw std::out_of_range("index out of range");
    }

    CV MatCSR::operator*(const CV& cv) const {
        assert(cv.m == n);
        CV cv2(m);
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = IA[i]; j < IA[i + 1]; j++) {
                cv2[i] += A[j] * cv[JA[j]];
            }
        }
        return cv2;
    }

    CV MatCSR::conjugateGradientSolve(CV& b, CV& x0, double maxError, uint64_t maxIterations) {
        CV x(x0.m);
        CV r(b.m);
        CV p(b.m);
        CV Ap(b.m);
        CV::copy(x, x0);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV::copy(p, r);
        double rTr = CV::squareSum(r);

        for (uint64_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
            double alpha = rTr / CV::dot(p, Ap);
            CV::addScaled(x, x, alpha, p);
            CV::addScaled(r, r, -alpha, Ap);
            double rTrNew = CV::squareSum(r);
            if (rTrNew < maxError * maxError) {
                printf("converged after %lu iterations\n", i);
                break;
            }
            double beta = rTrNew / rTr;
            CV::addScaled(p, r, beta, p);
            rTr = rTrNew;
        }
        return x;
    }

    MatCSR MatCSR::operator=(const MatCSR& other) {
        if (m != other.m || n != other.n || nnz != other.nnz) {
            free(IA);
            free(JA);
            free(A);
            m = other.m;
            n = other.n;
            nnz = other.nnz;
            IA = (uint64_t*) calloc(m + 1, sizeof(uint64_t));
            JA = (uint64_t*) calloc(nnz, sizeof(uint64_t));
            A = (double*) calloc(nnz, sizeof(double));
        }
        for (uint64_t i = 0; i < m + 1; i++) {
            IA[i] = other.IA[i];
        }
        for (uint64_t i = 0; i < nnz; i++) {
            JA[i] = other.JA[i];
            A[i] = other.A[i];
        }
        return *this;
    }
}