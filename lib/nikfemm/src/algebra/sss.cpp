#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <math.h>

#include <constants.hpp>

#include "sss.hpp"

#include "../utils/utils.hpp"

namespace nikfemm {
    MatSSS::MatSSS(MatCOO& coo) {
        m = coo.m;
        n = coo.n;
        nnz = coo.elems.size();

        row_ptr = new uint32_t[m + 1]();
        col_ind = new uint32_t[nnz - m];
        val = new double[nnz - m];
        diag = new double[m];

        uint32_t i = 0;
        for (auto const& [key, value] : coo.elems) {
            uint32_t _m = key >> 32;
            uint32_t _n = key & 0xFFFFFFFF;
            if (_m == _n) {
                diag[_m] = value;
            } else if (_m  > _n) {  // < for upper triangular, > for lower triangular
                col_ind[i] = _n;
                val[i] = value;
                row_ptr[_m + 1]++;
                i++;
            }
        }

        for (uint32_t i = 0; i < m; i++) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    MatSSS::~MatSSS() {
        delete[] row_ptr;
        delete[] col_ind;
        delete[] val;
        delete[] diag;
    }

    void MatSSS::printSSS() {
        printf("m: %lu, n: %lu, nnz: %lu\n", m, n, nnz);
        printf("row_ptr: ");
        for (uint32_t i = 0; i < m + 1; i++) {
            printf("%lu ", row_ptr[i]);
        }
        printf("\n");
        printf("col_ind: ");
        for (uint32_t i = 0; i < nnz - m; i++) {
            printf("%lu ", col_ind[i]);
        }
        printf("\n");
        printf("val: ");
        for (uint32_t i = 0; i < nnz - m; i++) {
            printf("%.1f ", val[i]);
        }
        printf("\n");
        printf("diag: ");
        for (uint32_t i = 0; i < m; i++) {
            printf("%.1f ", diag[i]);
        }
        printf("\n");
    }

    void MatSSS::print() {
        printf("m: %lu, n: %lu, nnz: %lu\n", m, n, nnz);
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < n; j++) {
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    void MatSSS::write_to_file(const char* filename) {
        FILE* fp = fopen(filename, "w");
        fprintf(fp, "%lu %lu %lu\n", m, n, nnz);
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < n; j++) {
                fprintf(fp, "%.17g ", (*this)(i, j));
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    double MatSSS::operator()(uint32_t i, uint32_t j) const {
        if (i == j) {
            return diag[i];
        } else if (i > j) {  // < for upper triangular, > for lower triangular
            uint32_t k = row_ptr[i];
            while (k < row_ptr[i + 1]) {
                if (col_ind[k] == j) {
                    return val[k];
                }
                k++;
            }
            return 0;
        } else {
            return 0;
        }
    }

    CV MatSSS::getInverseDiagonal() const {
        CV invDiag(m);
        for (uint32_t i = 0; i < m; i++) {
            invDiag[i] = 1 / diag[i];
        }
        return invDiag;
    }

    void MatSSS::preconditionedConjugateGradientSolve(CV& b, CV& x, double maxError, uint32_t maxIterations) {
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
}