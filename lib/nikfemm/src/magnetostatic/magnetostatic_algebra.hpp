#ifndef NIK_MAGNETOSTATIC_ALGEBRA_HPP
#define NIK_MAGNETOSTATIC_ALGEBRA_HPP

#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticNonLinearTerm {
        double linear_coefficient;
        const MagnetostaticProp* nonlinear_coefficient;
        uint32_t nonlinear_coefficient_element_index;
    };

    struct MagnetostaticNonLinearExpression {
        std::vector<MagnetostaticNonLinearTerm> terms;
    };

    struct MagnetostaticMatCSRSymmetric : virtual MatCSRSymmetric {
        MagnetostaticNonLinearExpression* diag_expr;  // these are used to update the stiffness matrix in an efficient way from the magnetic induction
        MagnetostaticNonLinearExpression* expr;

        MagnetostaticMatCSRSymmetric(MatCOO<MagnetostaticNonLinearExpression>& coo);
        ~MagnetostaticMatCSRSymmetric();

        void updateNonLinearCoefficients(std::vector<Vector>& B);
    };
}

#endif