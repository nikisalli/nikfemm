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
    MatSSS::MatSSS(const MatCOO& coo) {
        assert(coo.m == coo.n);
        m = coo.m;
        // count number of non-zero elements on coo diagonal
        uint32_t nnz_diag = 0;
        
        std::vector<std::pair<uint64_t, double>> elems;
        elems.reserve(coo.elems.size());

        for (auto const& [key, val] : coo.elems) {
            elems.push_back({key, val});
            nnz_diag += (key >> 32) == (key & 0xFFFFFFFF) * (val != 0);
        }
        printf("nnz_diag: %lu\n", nnz_diag);

        // sort elems by key
        std::sort(elems.begin(), elems.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        nnz = coo.elems.size() - nnz_diag;

        row_ptr = new uint32_t[m + 1]();
        col_ind = new uint32_t[nnz];
        val = new double[nnz];  // fix
        diag = new double[m];

        uint32_t i = 0;
        for (auto const& [key, value] : elems) {
            uint32_t _m = key >> 32;
            uint32_t _n = key & 0xFFFFFFFF;
            printf("key: %lu, _m: %lu, _n: %lu, value: %f\n", key, _m, _n, value);
            if (_m == _n) {
                diag[_m] = value;
            } else if (_m  < _n) {  // < for upper triangular, > for lower triangular
                col_ind[i] = _m;  // m and n are swapped
                val[i] = value;
                row_ptr[_n + 1]++;
                i++;
            }
        }

        for (uint32_t i = 0; i < m; i++) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    MatSSS::MatSSS(const MatSSS& sss) {
        copy(*this, sss);
    }

    MatSSS::~MatSSS() {
        delete[] row_ptr;
        delete[] col_ind;
        delete[] val;
        delete[] diag;
    }

    void MatSSS::printSSS() {
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
        printf("val: ");
        for (uint32_t i = 0; i < nnz; i++) {
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
        printf("m: %lu, nnz: %lu\n", m, nnz);
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                if (i <= j) {
                    printf("%.1f ", (*this)(i, j));
                } else {
                    printf("%.1f ", (*this)(j, i));
                }
            }
            printf("\n");
        }
    }

    void MatSSS::copy(MatSSS& result, const MatSSS& mat) {
        result.m = mat.m;
        result.nnz = mat.nnz;

        result.row_ptr = new uint32_t[result.m + 1];
        result.col_ind = new uint32_t[result.nnz];
        result.val = new double[result.nnz];
        result.diag = new double[result.m];

        memcpy(result.row_ptr, mat.row_ptr, (result.m + 1) * sizeof(uint32_t));
        memcpy(result.col_ind, mat.col_ind, (result.nnz) * sizeof(uint32_t));
        memcpy(result.val, mat.val, result.nnz * sizeof(double));
        memcpy(result.diag, mat.diag, result.m * sizeof(double));
    }

    void MatSSS::write_to_file(const char* filename) {
        FILE* fp = fopen(filename, "w");
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                if (i <= j) {
                    fprintf(fp, "%.17g ", (*this)(i, j));
                } else {
                    fprintf(fp, "%.17g ", (*this)(j, i));
                }
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    double MatSSS::operator()(const uint32_t i, const uint32_t j) const {
        if (i == j) {
            return diag[i];
        } else if (i < j) {  // < for upper triangular, > for lower triangular
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_ind[k] == j) {
                    return val[k];
                }
            }
            return 0;
        } else {
            return (this->operator()(j, i));
        }
    }

    void MatSSS::multSSORPreconditioner(CV& result, const CV& cv, double omega) {
        printf("multSSORPreconditioner\n");
        result.print();
        cv.print();
        double temp = omega * (2 - omega);
        for (uint32_t i = 0; i < m; i++) {
            result[i] = cv[i] * temp;
        }

        // invert lower triangle
        for (int64_t i = m - 1; i >= 0; i--) {
            result[i] /= diag[i];
            printf("i: %lu, diag[i]: %.1f, result[i]: %.1f\n", i, diag[i], result[i]);
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                printf("- k: %lu, col_ind[k]: %lu, val[k]: %.1f\n", k, col_ind[k], val[k]);
                result[col_ind[k]] -= val[k] * result[i] * omega;
            }
        }

        result.print();

        for (uint32_t i = 0; i < m; i++) {
            result[i] *= diag[i];
        }

        // invert upper triangle
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                result[i] -= val[k] * result[col_ind[k]] * omega;
            }
            result[i] /= diag[i];
        }
        result.print();
    }

    void MatSSS::preconditionedJacobiConjugateGradientSolver(CV& b, CV& x, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV P(r.m);
        for (uint32_t i = 0; i < m; i++) {
            P[i] = 1 / diag[i];
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
                printf("converged after %lu iterations\n", i);
#ifdef DEBUG_PRINT
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

    void MatSSS::preconditionedSSORConjugateGradientSolver(CV& b, CV& x, double omega, double maxError, uint32_t maxIterations) {
        CV r(b.m);
        CV::mult(r, *this, x);
        CV::sub(r, b, r);
        CV z(b.m);
        MatSSS::multSSORPreconditioner(z, r, omega);
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
            MatSSS::multSSORPreconditioner(z, r, omega);
            double beta = CV::dot(r, z) / rTzold;
            CV::addScaled(p, z, beta, p);
        }
    }
}