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
    BaseCSR::BaseCSR() {
        row_ptr = std::vector<uint32_t>();
        col_ind = std::vector<uint32_t>();
        val = std::vector<double>();
        diag = std::vector<double>();
        nnz = 0;
        m = 0;
    }

    void BaseCSR::printCSR() {
        printf("m: %u, nnz: %u\n", m, nnz);
        printf("row_ptr: ");
        for (uint32_t i = 0; i < m + 1; i++) {
            printf("%u ", row_ptr[i]);
        }
        printf("\n");
        printf("col_ind: ");
        for (uint32_t i = 0; i < nnz; i++) {
            printf("%u ", col_ind[i]);
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

    BaseCSR::BaseCSR(BuildMatCOO<double>& coo) {
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

        row_ptr = std::vector<uint32_t>(m + 1, 0);
        col_ind = std::vector<uint32_t>(nnz);
        val = std::vector<double>(nnz);
        diag = std::vector<double>(m);

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

        row_ptr = std::vector<uint32_t>(csr.row_ptr);
        col_ind = std::vector<uint32_t>(csr.col_ind);
        val = std::vector<double>(csr.val);
        diag = std::vector<double>(csr.diag);
    }

    BaseCSR::~BaseCSR() {
        
    }

    // MatCSRSymmetric

    MatCSRSymmetric::MatCSRSymmetric() : BaseCSR() {}

    void MatCSRSymmetric::print() {
        // iterate over CSR elements
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
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = 0; j < m; j++) {
                fprintf(f, "%.17g ", (*this)(i, j));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    // MatCSRLowerTri

    MatCSRLowerTri::MatCSRLowerTri() : BaseCSR() {}

    void MatCSRLowerTri::print() {
        // iterate over CSR elements
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
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = 0; j < m; j++) {
                fprintf(f, "%.17g ", (*this)(i, j));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    // MatCSRUpperTri

    MatCSRUpperTri::MatCSRUpperTri() : BaseCSR() {}

    void MatCSRUpperTri::print() {
        // iterate over CSR elements
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
        for (uint64_t i = 0; i < m; i++) {
            for (uint64_t j = 0; j < m; j++) {
                fprintf(f, "%.17g ", (*this)(i, j));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
}