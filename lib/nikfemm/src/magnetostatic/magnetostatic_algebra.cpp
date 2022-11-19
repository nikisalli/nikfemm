#include <algorithm>

#include "magnetostatic_algebra.hpp"

namespace nikfemm {
    MagnetostaticMatCSRSymmetric::MagnetostaticMatCSRSymmetric(MatCOO<MagnetostaticNonLinearExpression>& coo) {
        m = coo.m;

        std::vector<std::pair<uint64_t, MagnetostaticNonLinearExpression>> elems;
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
        expr = new MagnetostaticNonLinearExpression[nnz];
        diag_expr = new MagnetostaticNonLinearExpression[m];

        uint32_t i = 0;
        for (auto const& [key, value] : elems) {
            uint32_t _m = key >> 32;
            uint32_t _n = key & 0xFFFFFFFF;
            if (_m == _n) {
                // diag[_m] = value;
                diag_expr[_m] = value;
            } else {
                col_ind[i] = _n;
                // val[i] = value;
                expr[i] = value;
                row_ptr[_m + 1]++;
                i++;
            }
        }

        for (uint32_t i = 0; i < m; i++) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    MagnetostaticMatCSRSymmetric::~MagnetostaticMatCSRSymmetric() {
        delete[] expr;
        delete[] diag_expr;
    }

    void MagnetostaticMatCSRSymmetric::updateMu(std::vector<const MagnetostaticProp*>& props, std::vector<float>& mu, std::vector<Vector>& B) {
        assert(mu.size() == B.size());
        for (uint32_t i = 0; i < B.size(); i++) {
            float Bmag = B[i].magnitude();
            mu[i] += (props[i]->getMu(Bmag) - mu[i]) * 0.1;
        }
    }

    void MagnetostaticMatCSRSymmetric::updateMat(std::vector<float>& mu) {
        // diagonal
        for (uint32_t i = 0; i < m; i++) {
            diag[i] = 0;
            for (auto term : diag_expr[i].terms) {
                if (term.is_boundary_condition) {
                    diag[i] = term.linear_coefficient;
                } else {
                    diag[i] += term.linear_coefficient / mu[term.nonlinear_coefficient_element_index];
                    // printf("diag[%d] += %.17g / mu[%d] (%.17g) = %.17g\n", i, term.linear_coefficient, term.nonlinear_coefficient_element_index, mu[term.nonlinear_coefficient_element_index], diag[i]);
                }
            }
        }
        // off-diagonal
        for (uint32_t i = 0; i < nnz; i++) {
            val[i] = 0;
            for (auto term : expr[i].terms) {
                if (term.is_boundary_condition) {
                    val[i] = term.linear_coefficient;
                } else {
                    val[i] += term.linear_coefficient / mu[term.nonlinear_coefficient_element_index];
                    // printf("val[%d] += %.17g / mu[%d] (%.17g) = %.17g\n", i, term.linear_coefficient, term.nonlinear_coefficient_element_index, mu[term.nonlinear_coefficient_element_index], diag[i]);
                }
            }
        }
    }
}