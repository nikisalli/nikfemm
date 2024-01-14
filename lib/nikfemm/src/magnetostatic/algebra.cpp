#include <algorithm>

#include "algebra.hpp"

namespace nikfemm {
    double MagnetostaticNonLinearExpression::evaluate(std::vector<float>& mu) {
        double result = constant;
        for (auto const& term : terms) {
            result += term.linear_coefficient / mu[term.nonlinear_coefficient_element_index];
        }
        return result;
    }

    void MagnetostaticNonLinearExpression::setToConstant(double constant) {
        this->constant = constant;
        terms.clear();
    }

    void MagnetostaticNonLinearExpression::addToConstant(double constant) {
        nassert(isLinear(), "Cannot add to constant of non-linear expression");
        this->constant += constant;
    }

    void MagnetostaticNonLinearExpression::addTerm(double linear_coefficient, uint32_t nonlinear_coefficient_element_index) {
        terms.push_back(MagnetostaticNonLinearTerm{linear_coefficient, nonlinear_coefficient_element_index});
    }

    MagnetostaticMatCSRSymmetric::MagnetostaticMatCSRSymmetric(BuildMatCOO<MagnetostaticNonLinearExpression>& coo) {
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

        row_ptr = std::vector<uint32_t>(m + 1, 0);
        col_ind = std::vector<uint32_t>(nnz);
        val = std::vector<double>(nnz);
        diag = std::vector<double>(m);
        expr = std::vector<MagnetostaticNonLinearExpression>(nnz);
        diag_expr = std::vector<MagnetostaticNonLinearExpression>(m);

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

    }

    void MagnetostaticMatCSRSymmetric::updateFromMu(std::vector<float>& mu) {
        // diagonal
        for (uint32_t i = 0; i < m; i++) {
            double old_val = diag[i];
            diag[i] = diag_expr[i].evaluate(mu);
        }
        // off-diagonal
        for (uint32_t i = 0; i < nnz; i++) {
            val[i] = expr[i].evaluate(mu);
        }
    }

    MagnetostaticCV::MagnetostaticCV(uint32_t size) {
        val = std::vector<double>(size);
        expr = std::vector<MagnetostaticNonLinearExpression>(size);
    }

    MagnetostaticCV::~MagnetostaticCV() {

    }

    void MagnetostaticCV::updateFromMu(std::vector<float>& mu) {
        for (uint32_t i = 0; i < val.size(); i++) {
            val[i] = expr[i].evaluate(mu);
        }
    }
}