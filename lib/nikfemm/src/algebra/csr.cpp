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
    void BaseCSR::printCSR() {
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
        printf("diag: ");
        for (uint32_t i = 0; i < m; i++) {
            printf("%.1f ", diag[i]);
        }
        printf("\n");
    }

    BaseCSR::BaseCSR(MatCOO& coo) {
        m = coo.m;

        std::vector<std::pair<uint64_t, double>> elems;
        elems.reserve(coo.elems.size());

        for (auto const& [key, value] : coo.elems) {
            elems.push_back(std::make_pair(key, value));
        }

        // sort elems by key
        std::sort(elems.begin(), elems.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        nnz = elems.size() - m;

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
            } else {
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

    BaseCSR::BaseCSR(const BaseCSR& csr) {
        m = csr.m;
        nnz = csr.nnz;

        row_ptr = new uint32_t[m + 1];
        col_ind = new uint32_t[nnz];
        val = new double[nnz];
        diag = new double[m];

        memcpy(row_ptr, csr.row_ptr, (m + 1) * sizeof(uint32_t));
        memcpy(col_ind, csr.col_ind, nnz * sizeof(uint32_t));
        memcpy(val, csr.val, nnz * sizeof(double));
        memcpy(diag, csr.diag, m * sizeof(double));
    }

    BaseCSR::~BaseCSR() {
        delete[] row_ptr;
        delete[] col_ind;
        delete[] val;
        delete[] diag;
    }

    // MatCSRSymmetric

    void MatCSRSymmetric::print() {
        // iterate over CSR elements
        uint32_t idx = 0;
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                printf("%.1f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    double MatCSRSymmetric::operator()(uint32_t i, uint32_t j) const {
        if (i == j) {
            return diag[i];
        } else if (i < j) {
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_ind[k] == j) {
                    return val[k];
                }
            }
        } else {
            return (*this)(j, i);
        }
        return 0.0;
    }

    void MatCSRSymmetric::write_to_file(const char *filename) {
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

    // MatCSRLowerTri

    void MatCSRLowerTri::print() {
        // iterate over CSR elements
        uint32_t idx = 0;
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                printf("%.4f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    double MatCSRLowerTri::operator()(uint32_t i, uint32_t j) const {
        if (i == j) {
            return diag[i];
        } else {
            for (uint32_t k = row_ptr[j]; k < row_ptr[j + 1]; k++) {
                if (col_ind[k] == i) {
                    return val[k];
                }
            }
        }
        return 0.0;
    }

    void MatCSRLowerTri::write_to_file(const char *filename) {
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

    // MatCSRUpperTri

    void MatCSRUpperTri::print() {
        // iterate over CSR elements
        uint32_t idx = 0;
        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < m; j++) {
                printf("%.4f ", (*this)(i, j));
            }
            printf("\n");
        }
    }

    double MatCSRUpperTri::operator()(uint32_t i, uint32_t j) const {
        if (i == j) {
            return diag[i];
        } else {
            for (uint32_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_ind[k] == j) {
                    return val[k];
                }
            }
        }
        return 0.0;
    }

    void MatCSRUpperTri::write_to_file(const char *filename) {
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
}