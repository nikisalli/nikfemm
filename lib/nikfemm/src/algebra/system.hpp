#ifndef NIK_SYSTEM_HPP
#define NIK_SYSTEM_HPP

#include "../algebra/csr.hpp"
#include "../algebra/coo.hpp"

namespace nikfemm {
    template <typename WeightType>
    struct System {
        MatCOOSymmetric<WeightType> A;
        std::vector<WeightType> b;

        // constructor
        System(uint32_t m) : A(m), b(m) {}
        System(MatCOOSymmetric<WeightType> A, std::vector<double> b) : A(A), b(b) {}
        System(System<WeightType> const& other) : A(other.A), b(other.b) {}

        void addDirichletBoundaryCondition(uint32_t id, WeightType value) {
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
                    b[n] -= elem.second * value;
                    elem.second = 0;
                }
                if (n == id) { // column index equal to id
                    b[m] -= elem.second * value;
                    elem.second = 0;
                }
            }

            // coo.elems[MatCOOSymmetric<int>::get_key(id, id)].setToConstant(1);
            A(id, id) = 1;
            b[id] = value;
        }
    };
}

#endif