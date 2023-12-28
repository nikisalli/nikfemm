#include <nikfemm.hpp>

#define BOOST_TEST_MODULE VsidCommonTest
#include <boost/test/unit_test.hpp>

using namespace nikfemm;

BOOST_AUTO_TEST_SUITE( csr_tests )

BOOST_AUTO_TEST_CASE( csr_index )
{
    BuildMatCOO<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    MatCSRSymmetric csr(mat);

    BOOST_CHECK_EQUAL(csr(0, 0), 1);
    BOOST_CHECK_EQUAL(csr(0, 1), 2);
    BOOST_CHECK_EQUAL(csr(0, 2), 3);
    BOOST_CHECK_EQUAL(csr(1, 1), 4);
    BOOST_CHECK_EQUAL(csr(1, 2), 5);
    BOOST_CHECK_EQUAL(csr(2, 2), 6);
}

BOOST_AUTO_TEST_CASE( csr_cv_mult )
{
    BuildMatCOO<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    MatCSRSymmetric csr(mat);

    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);

    CV::mult(y, csr, x);

    // wolfram alpha: {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}} * {1, 2, 3}

    BOOST_CHECK_EQUAL(y[0], 14);
    BOOST_CHECK_EQUAL(y[1], 25);
    BOOST_CHECK_EQUAL(y[2], 31);
}

BOOST_AUTO_TEST_CASE( csr_conjugate_gradient_solve )
{
    BuildMatCOO<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    MatCSRSymmetric csr(mat);

    CV x(3);
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;

    CV b(3);
    b[0] = 14;
    b[1] = 32;
    b[2] = 54;

    preconditionedSSORConjugateGradientSolver(csr, b, x, 1.5, 1e-6, 1000);

    // wolfram alpha: solve {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}} * {x, y, z} = {14, 32, 54}
    // solution: {x, y, z} = {26, 0, -4}

    // check result inside tolerance
    BOOST_CHECK_SMALL(x[0] - 26, 1e-6);
    BOOST_CHECK_SMALL(x[1] - 0, 1e-6);
    BOOST_CHECK_SMALL(x[2] - -4, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( cv_tests )

BOOST_AUTO_TEST_CASE( cv_index )
{
    CV cv(3);
    cv[0] = 1;
    cv[1] = 2;
    cv[2] = 3;

    BOOST_CHECK_EQUAL(cv[0], 1);
    BOOST_CHECK_EQUAL(cv[1], 2);
    BOOST_CHECK_EQUAL(cv[2], 3);
}

BOOST_AUTO_TEST_CASE( cv_mult )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    CV z(3);
    
    CV::mult(z, x, y);

    BOOST_CHECK_EQUAL(z[0], 4);
    BOOST_CHECK_EQUAL(z[1], 10);
    BOOST_CHECK_EQUAL(z[2], 18);
}

BOOST_AUTO_TEST_CASE( cv_mult_scalar )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    
    CV::mult(y, 2, x);

    BOOST_CHECK_EQUAL(y[0], 2);
    BOOST_CHECK_EQUAL(y[1], 4);
    BOOST_CHECK_EQUAL(y[2], 6);
}

BOOST_AUTO_TEST_CASE( cv_add )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    CV z(3);
    
    CV::add(z, x, y);

    BOOST_CHECK_EQUAL(z[0], 5);
    BOOST_CHECK_EQUAL(z[1], 7);
    BOOST_CHECK_EQUAL(z[2], 9);
}

BOOST_AUTO_TEST_CASE( cv_sub )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    CV z(3);
    
    CV::sub(z, x, y);

    BOOST_CHECK_EQUAL(z[0], -3);
    BOOST_CHECK_EQUAL(z[1], -3);
    BOOST_CHECK_EQUAL(z[2], -3);
}

BOOST_AUTO_TEST_CASE( cv_div )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    CV z(3);
    
    CV::div(z, x, y);

    BOOST_CHECK_EQUAL(z[0], 0.25);
    BOOST_CHECK_EQUAL(z[1], 0.4);
    BOOST_CHECK_EQUAL(z[2], 0.5);
}

BOOST_AUTO_TEST_CASE( cv_div_scalar )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    
    CV::div(y, x, 2);

    BOOST_CHECK_EQUAL(y[0], 0.5);
    BOOST_CHECK_EQUAL(y[1], 1);
    BOOST_CHECK_EQUAL(y[2], 1.5);
}

BOOST_AUTO_TEST_CASE( cv_sub_scalar )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    
    CV::sub(y, x, 2);

    BOOST_CHECK_EQUAL(y[0], -1);
    BOOST_CHECK_EQUAL(y[1], 0);
    BOOST_CHECK_EQUAL(y[2], 1);
}

BOOST_AUTO_TEST_CASE( cv_add_scalar )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    
    CV::add(y, x, 2);

    BOOST_CHECK_EQUAL(y[0], 3);
    BOOST_CHECK_EQUAL(y[1], 4);
    BOOST_CHECK_EQUAL(y[2], 5);
}

BOOST_AUTO_TEST_CASE( cv_dot )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    float z = CV::dot(x, y);

    BOOST_CHECK_EQUAL(z, 32);
}

BOOST_AUTO_TEST_CASE( cv_add_scaled )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    CV z(3);
    
    CV::addScaled(z, x, 2, y);

    BOOST_CHECK_EQUAL(z[0], 9);
    BOOST_CHECK_EQUAL(z[1], 12);
    BOOST_CHECK_EQUAL(z[2], 15);
}

BOOST_AUTO_TEST_CASE( cv_square_sum )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    float z = CV::squareSum(x);

    BOOST_CHECK_EQUAL(z, 14);
}

BOOST_AUTO_TEST_CASE( cv_copy )
{
    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    
    CV::copy(y, x);

    BOOST_CHECK_EQUAL(y[0], 1);
    BOOST_CHECK_EQUAL(y[1], 2);
    BOOST_CHECK_EQUAL(y[2], 3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( current_density_tests )

BOOST_AUTO_TEST_CASE( current_density_dirichlet_boundary_conditions )
{
    BuildMatCOO<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    CV b(3);
    b[0] = 14;
    b[1] = 32;
    b[2] = 54;

    CurrentDensitySystem system = {
        mat,
        b
    };

    printf("A without Dirichlet boundary conditions:\n");
    MatCSRSymmetric csr_1(system.A);
    csr_1.print();
    printf("b without Dirichlet boundary conditions:\n");
    system.b.print();

    system.addDirichletBoundaryCondition(0, 1.0);
    system.addDirichletBoundaryCondition(1, 2.0);

    MatCSRSymmetric csr(system.A);

    printf("A with Dirichlet boundary conditions:\n");
    csr.print();
    printf("b with Dirichlet boundary conditions:\n");
    system.b.print();

    CV x(3);
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;

    preconditionedSSORConjugateGradientSolver(csr, system.b, x, 1.5, 1e-6, 1000);

    // wolfram alpha: solve {{1, 0, 0}, {0, 4, 5}, {0, 5, 6}} * {x, y, z} = {1, 30, 51}
    // solution: {x, y, z} = {1, 75, -54}

    // check result inside tolerance
    BOOST_CHECK_SMALL(x[0] - 1, 1e-6);
    BOOST_CHECK_SMALL(x[1] - 75, 1e-6);
    BOOST_CHECK_SMALL(x[2] - -54, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()