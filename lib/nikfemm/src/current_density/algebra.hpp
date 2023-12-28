#ifndef NIK_CURRENT_DENSITY_ALGEBRA_HPP
#define NIK_CURRENT_DENSITY_ALGEBRA_HPP

#include "../algebra/csr.hpp"
#include "../algebra/simple_vector.hpp"
#include "../algebra/build_coo.hpp"
#include "properties.hpp"

namespace nikfemm {
    struct CurrentDensitySystem {
        BuildMatCOO<double> A;
        CV b;

        void addDirichletBoundaryCondition(uint32_t id, double value) {
            // https://community.freefem.org/t/implementation-of-dirichlet-boundary-condition-when-tgv-1/113
            // this function lets you set a Dirichlet boundary condition on a node

            // every element that has a row or column index equal to id is set to zero
            // for every element that has a column index equal to id, its value multiplied 
            // by the value of the boundary condition is subtracted from the corresponding element in the b vector
            // care must be taken because the matrix is symmetric but it is stored as upper triangular to save memory

            for (auto& elem : A.elems) {
                uint32_t m = elem.first >> 32;
                uint32_t n = elem.first & 0xFFFFFFFF;

                // skip diagonal elements
                if (m == n) continue;

                if (m == id) { // row index equal to id
                    b.val[n] -= elem.second * value;
                    elem.second = 0;
                }
                if (n == id) { // column index equal to id
                    b.val[m] -= elem.second * value;
                    elem.second = 0;
                }
            }

            // coo.elems[BuildMatCOO<int>::get_key(id, id)].setToConstant(1);
            A(id, id) = 1;
            b.val[id] = value;
        }
    };
}

#endif