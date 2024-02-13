#ifndef NIK_MAGNETOSTATIC_ALGEBRA_HPP
#define NIK_MAGNETOSTATIC_ALGEBRA_HPP

#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"

namespace nikfemm {
    struct NonLinearTerm {
        double linear_coefficient;
        uint32_t nonlinear_coefficient_element_index;
    };

    struct NonLinearExpression {
        private:
        double constant = 0;
        std::vector<NonLinearTerm> terms;

        public:
        inline bool isLinear() {
            return terms.size() == 0;
        }
        double evaluate(std::vector<double>& mu);

        NonLinearExpression& operator=(const double& other) {
            constant = other;
            terms.clear();
            return *this;
        }

        NonLinearExpression& operator=(const NonLinearExpression& other) {
            constant = other.constant;
            terms = other.terms;
            return *this;
        }

        NonLinearExpression& operator+=(const double& other) {
            constant += other;
            return *this;
        }

        NonLinearExpression& operator+=(const NonLinearTerm& other) {
            terms.push_back(other);
            return *this;
        }

        NonLinearExpression& operator+=(const NonLinearExpression& other) {
            constant += other.constant;
            for (auto const& term : other.terms) {
                terms.push_back(term);
            }
            return *this;
        }

        NonLinearExpression& operator-=(const double& other) {
            constant -= other;
            return *this;
        }

        NonLinearExpression& operator-=(const NonLinearTerm& other) {
            terms.push_back(NonLinearTerm{-other.linear_coefficient, other.nonlinear_coefficient_element_index});
            return *this;
        }

        NonLinearExpression& operator-=(const NonLinearExpression& other) {
            constant -= other.constant;
            for (auto const& term : other.terms) {
                terms.push_back(NonLinearTerm{-term.linear_coefficient, term.nonlinear_coefficient_element_index});
            }
            return *this;
        }

        NonLinearExpression operator*(const double& other) {
            NonLinearExpression result;
            result.constant = constant * other;
            for (auto const& term : terms) {
                result.terms.push_back(NonLinearTerm{term.linear_coefficient * other, term.nonlinear_coefficient_element_index});
            }
            return result;
        }
    };

    struct NonLinearMatCSRSymmetric : virtual MatCSRSymmetric {
        std::vector<NonLinearExpression*> diag_expr;  // these are used to update the stiffness matrix in an efficient way from the magnetic induction
        std::vector<NonLinearExpression*> expr;
        // NonLinearExpression* diag_expr;  // these are used to update the stiffness matrix in an efficient way from the magnetic induction
        // NonLinearExpression* expr;

        NonLinearMatCSRSymmetric(MatCOOSymmetric<NonLinearExpression>& coo);

        void evaluate(std::vector<double>& mu);
    };
}

#endif