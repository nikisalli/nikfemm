#include <nikfemm.hpp>
#include <stdio.h>
#include <stdlib.h>

using namespace nikfemm;

int main(int argc, char** argv) {
    MatCOO coo(8, 8);
    //   182     0    36     0     0    99    80     0
    //     0     4     0     0     0     0     0     0
    //    36     0   211     0     0    99    90     0
    //     0     0     0    16     0     0     0     0
    //     0     0     0     0   169     0     0   156
    //    99     0    99     0     0   157     0     0
    //    80     0    90     0     0     0   149     0
    //     0     0     0     0   156     0     0   208
    coo.set_elem(0, 0, 182);
    coo.set_elem(0, 2, 36);
    coo.set_elem(0, 5, 99);
    coo.set_elem(0, 6, 80);
    coo.set_elem(1, 1, 4);
    coo.set_elem(2, 0, 36);
    coo.set_elem(2, 2, 211);
    coo.set_elem(2, 5, 99);
    coo.set_elem(2, 6, 90);
    coo.set_elem(3, 3, 16);
    coo.set_elem(4, 4, 169);
    coo.set_elem(4, 7, 156);
    coo.set_elem(5, 0, 99);
    coo.set_elem(5, 2, 99);
    coo.set_elem(5, 5, 157);
    coo.set_elem(6, 0, 80);
    coo.set_elem(6, 2, 90);
    coo.set_elem(6, 6, 149);
    coo.set_elem(7, 4, 156);
    coo.set_elem(7, 7, 208);

    CV b(8);
    b.set_elem(0, 1);
    b.set_elem(1, 2);
    b.set_elem(2, 3);
    b.set_elem(3, 4);
    b.set_elem(4, 5);
    b.set_elem(5, 6);
    b.set_elem(6, 7);
    b.set_elem(7, 8);

    CV x(8);
    
    MatCSR csr(coo);
    csr.preconditionedJacobiConjugateGradientSolver(x, b, 1e-6, 10);
    // csr.preconditionedSSORConjugateGradientSolver(b, x, 1, 1e-7, 100);

    csr.print();
}