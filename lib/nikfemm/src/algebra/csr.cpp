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
        IA = (uint64_t*) calloc(m + 1, sizeof(uint64_t));
        JA = (uint64_t*) calloc(nnz, sizeof(uint64_t));
        A = (double*) calloc(nnz, sizeof(double));

        IA[0] = 0;

        struct ElemCoo {
            uint64_t m;
            uint64_t n;
            double val;
        };

        std::vector<ElemCoo> elems_coo;

        for (auto elem : coo.elems) {
            elems_coo.push_back(
                ElemCoo {
                    elem.first >> 32,
                    elem.first & 0xFFFFFFFF,
                    elem.second
                }
            );
        }

        // sort by m and n
        std::sort(elems_coo.begin(), elems_coo.end(), [](ElemCoo& e1, ElemCoo& e2) {
            if (e1.m == e2.m) {
                return e1.n < e2.n;
            } else {
                return e1.m < e2.m;
            }
        });

        for (uint64_t i = 0; i < nnz; i++) {
            JA[i] = elems_coo[i].n;
            A[i] = elems_coo[i].val;
            IA[elems_coo[i].m + 1]++;
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
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    double MatCSR::operator()(uint64_t i, uint64_t j) const {
        // bounds check
        if (i >= m || j >= n) {
            nexit("MatCSR::operator(): index out of range");
        }
        uint64_t row_start = IA[i];
        uint64_t row_end = IA[i + 1];
        for (uint64_t k = row_start; k < row_end; k++) {
            if (JA[k] == j) {
                return A[k];
            }
        }
        return 0;
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

    void MatCSR::conjugateGradientSolve(CV& b, CV& x, double maxError, uint64_t maxIterations) {
        CV r(b.m);
        CV p(b.m);
        CV Ap(b.m);
        // printf("x0:\n");
        // x.print();
        CV::mult(r, *this, x);
        // printf("r:\n");
        // r.print();
        CV::sub(r, b, r);
        // printf("r:\n");
        // r.print();
        CV::copy(p, r);
        // printf("p:\n");
        // p.print();
        double rTr = CV::squareSum(r);
        // printf("rTr: %.1f\n", rTr);

        for (uint64_t i = 0; i < maxIterations; i++) {
            CV::mult(Ap, *this, p);
            // printf("Ap:\n");
            // Ap.print();
            // printf("p.Ap = %f\n", CV::dot(p, Ap));
            double pAp = CV::dot(p, Ap);
            if (fabs(pAp) < std::numeric_limits<double>::epsilon()) {
                pAp = std::numeric_limits<double>::epsilon();
                printf("warning: pAp is zero. approximating with epsilon");
            }
            double alpha = rTr / pAp;
            // printf("alpha: %.1f\n", alpha);
            CV::addScaled(x, x, alpha, p);
            // printf("x:\n");
            // x.print();
            CV::addScaled(r, r, -alpha, Ap);
            // printf("r:\n");
            // r.print();
            double rTrNew = CV::squareSum(r);
            // printf("rTrNew: %.1f\n", rTrNew);
            printf("iteration %lu, error: %f\n", i, sqrt(rTrNew));
            if (rTrNew < maxError * maxError) {
                printf("converged after %lu iterations\n", i);
                printf("x:\n");
                x.print();
                break;
            }
            double beta = rTrNew / rTr;
            // printf("beta: %.1f\n", beta);
            CV::addScaled(p, r, beta, p);
            // printf("p:\n");
            // p.print();
            rTr = rTrNew;
            // printf("rTr: %.1f\n", rTr);
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