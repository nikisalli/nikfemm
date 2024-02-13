#include <algorithm>

#include "non_linear.hpp"

namespace nikfemm {
    double NonLinearExpression::evaluate(std::vector<double>& mu) {
        double result = constant;
        for (auto const& term : terms) {
            result += term.linear_coefficient / mu[term.nonlinear_coefficient_element_index];
        }
        return result;
    }

    NonLinearMatCSRSymmetric::NonLinearMatCSRSymmetric(MatCOOSymmetric<NonLinearExpression>& coo) {
        m = coo.m;

        std::vector<uint64_t> keys(coo.elems.size());

        std::transform(coo.elems.begin(), coo.elems.end(), keys.begin(), [](const auto& elem) {
            return elem.first;
        });

        // sort elems by key
        std::sort(keys.begin(), keys.end());

        nnz = coo.elems.size() - m;

        row_ptr = std::vector<uint32_t>(m + 1, 0);
        col_ind = std::vector<uint32_t>(nnz);
        val = std::vector<double>(nnz);
        diag = std::vector<double>(m);
        expr = std::vector<NonLinearExpression*>(nnz);
        diag_expr = std::vector<NonLinearExpression*>(m);

        uint32_t i = 0;
        for (auto const key : keys) {
            uint32_t _m = key >> 32;
            uint32_t _n = key & 0xFFFFFFFF;
            if (_m == _n) {
                // diag[_m] = value;
                // move the expression to diag_expr
                diag_expr[_m] = &coo.elems[key];
            } else {
                col_ind[i] = _n;
                // val[i] = value;
                expr[i] = &coo.elems[key];
                row_ptr[_m + 1]++;
                i++;
            }
        }

        for (uint32_t i = 0; i < m; i++) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    void NonLinearMatCSRSymmetric::evaluate(std::vector<double>& mu) {
        // diagonal
        for (uint32_t i = 0; i < m; i++) {
            double old_val = diag[i];
            diag[i] = diag_expr[i]->evaluate(mu);
        }
        // off-diagonal
        for (uint32_t i = 0; i < nnz; i++) {
            val[i] = expr[i]->evaluate(mu);
        }
    }
}