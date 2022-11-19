#ifndef NIK_MAGNETOSTATIC_ALGEBRA_HPP
#define NIK_MAGNETOSTATIC_ALGEBRA_HPP

#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct MagnetostaticNonLinearTerm {
        double linear_coefficient;
        uint32_t nonlinear_coefficient_element_index;
        bool is_boundary_condition;
    };

    struct MagnetostaticNonLinearExpression {
        std::vector<MagnetostaticNonLinearTerm> terms;
    };

    struct MagnetostaticMatCSRSymmetric : virtual MatCSRSymmetric {
        MagnetostaticNonLinearExpression* diag_expr;  // these are used to update the stiffness matrix in an efficient way from the magnetic induction
        MagnetostaticNonLinearExpression* expr;

        MagnetostaticMatCSRSymmetric(MatCOO<MagnetostaticNonLinearExpression>& coo);
        ~MagnetostaticMatCSRSymmetric();

        void updateMu(std::vector<const MagnetostaticProp*>& props, std::vector<float>& mu, std::vector<Vector>& B);
        void updateMat(std::vector<float>& mu);
    };
}

#endif