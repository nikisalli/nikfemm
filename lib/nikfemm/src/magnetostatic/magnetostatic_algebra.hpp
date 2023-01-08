#ifndef NIK_MAGNETOSTATIC_ALGEBRA_HPP
#define NIK_MAGNETOSTATIC_ALGEBRA_HPP

#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/build_coo.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticNonLinearTerm {
        double linear_coefficient;
        uint32_t nonlinear_coefficient_element_index;
    };

    struct MagnetostaticNonLinearExpression {
        private:
        double constant = 0;
        std::vector<MagnetostaticNonLinearTerm> terms;

        public:
        inline bool isLinear() {
            return terms.size() == 0;
        }
        double evaluate(std::vector<float>& mu);
        void setToConstant(double constant);
        void setDirichlet();
        void unsetDirichlet();
        void addToConstant(double constant);
        void addTerm(double linear_coefficient, uint32_t nonlinear_coefficient_element_index);

        MagnetostaticNonLinearExpression& operator=(const MagnetostaticNonLinearExpression& other) {
            constant = other.constant;
            terms = other.terms;
            return *this;
        }

        MagnetostaticNonLinearExpression& operator+=(const MagnetostaticNonLinearExpression& other) {
            constant += other.constant;
            for (auto const& term : other.terms) {
                terms.push_back(term);
            }
            return *this;
        }

        MagnetostaticNonLinearExpression& operator-=(const MagnetostaticNonLinearExpression& other) {
            constant -= other.constant;
            for (auto const& term : other.terms) {
                terms.push_back(MagnetostaticNonLinearTerm{-term.linear_coefficient, term.nonlinear_coefficient_element_index});
            }
            return *this;
        }

        MagnetostaticNonLinearExpression operator*(const double& other) {
            MagnetostaticNonLinearExpression result;
            result.constant = constant * other;
            for (auto const& term : terms) {
                result.terms.push_back(MagnetostaticNonLinearTerm{term.linear_coefficient * other, term.nonlinear_coefficient_element_index});
            }
            return result;
        }
    };

    struct MagnetostaticMatCSRSymmetric : virtual MatCSRSymmetric {
        std::vector<MagnetostaticNonLinearExpression> diag_expr;  // these are used to update the stiffness matrix in an efficient way from the magnetic induction
        std::vector<MagnetostaticNonLinearExpression> expr;
        // MagnetostaticNonLinearExpression* diag_expr;  // these are used to update the stiffness matrix in an efficient way from the magnetic induction
        // MagnetostaticNonLinearExpression* expr;

        MagnetostaticMatCSRSymmetric(BuildMatCOO<MagnetostaticNonLinearExpression>& coo);
        ~MagnetostaticMatCSRSymmetric();

        void updateFromMu(std::vector<float>& mu);
    };

    struct MagnetostaticCV : virtual CV {
        std::vector<MagnetostaticNonLinearExpression> expr;

        MagnetostaticCV(uint32_t size);
        ~MagnetostaticCV();

        void updateFromMu(std::vector<float>& mu);
    };

    struct MagnetostaticSystem {
        BuildMatCOO<MagnetostaticNonLinearExpression> coo;
        MagnetostaticCV b;
    };
}

#endif