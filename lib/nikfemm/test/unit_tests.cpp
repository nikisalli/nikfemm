#include <nikfemm.hpp>

#define BOOST_TEST_MODULE VsidCommonTest
#include <boost/test/unit_test.hpp>

using namespace nikfemm;

BOOST_AUTO_TEST_SUITE( csr_tests )

BOOST_AUTO_TEST_CASE( csr_index )
{
    MatCOOSymmetric<double> mat(3);

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
    MatCOOSymmetric<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    MatCSRSymmetric csr(mat);

    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);

    mult(y, csr, x);

    // wolfram alpha: {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}} * {1, 2, 3}

    BOOST_CHECK_EQUAL(y[0], 14);
    BOOST_CHECK_EQUAL(y[1], 25);
    BOOST_CHECK_EQUAL(y[2], 31);
}

BOOST_AUTO_TEST_CASE( csr_conjugate_gradient_solve )
{
    MatCOOSymmetric<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    MatCSRSymmetric csr(mat);

    std::vector x(3);
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;

    std::vector b(3);
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
    std::vector cv(3);
    cv[0] = 1;
    cv[1] = 2;
    cv[2] = 3;

    BOOST_CHECK_EQUAL(cv[0], 1);
    BOOST_CHECK_EQUAL(cv[1], 2);
    BOOST_CHECK_EQUAL(cv[2], 3);
}

BOOST_AUTO_TEST_CASE( cv_mult )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    std::vector z(3);
    
    mult(z, x, y);

    BOOST_CHECK_EQUAL(z[0], 4);
    BOOST_CHECK_EQUAL(z[1], 10);
    BOOST_CHECK_EQUAL(z[2], 18);
}

BOOST_AUTO_TEST_CASE( cv_mult_scalar )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    
    mult(y, 2, x);

    BOOST_CHECK_EQUAL(y[0], 2);
    BOOST_CHECK_EQUAL(y[1], 4);
    BOOST_CHECK_EQUAL(y[2], 6);
}

BOOST_AUTO_TEST_CASE( cv_add )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    std::vector z(3);
    
    add(z, x, y);

    BOOST_CHECK_EQUAL(z[0], 5);
    BOOST_CHECK_EQUAL(z[1], 7);
    BOOST_CHECK_EQUAL(z[2], 9);
}

BOOST_AUTO_TEST_CASE( cv_sub )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    std::vector z(3);
    
    sub(z, x, y);

    BOOST_CHECK_EQUAL(z[0], -3);
    BOOST_CHECK_EQUAL(z[1], -3);
    BOOST_CHECK_EQUAL(z[2], -3);
}

BOOST_AUTO_TEST_CASE( cv_div )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    std::vector z(3);
    
    div(z, x, y);

    BOOST_CHECK_EQUAL(z[0], 0.25);
    BOOST_CHECK_EQUAL(z[1], 0.4);
    BOOST_CHECK_EQUAL(z[2], 0.5);
}

BOOST_AUTO_TEST_CASE( cv_div_scalar )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    
    div(y, x, 2);

    BOOST_CHECK_EQUAL(y[0], 0.5);
    BOOST_CHECK_EQUAL(y[1], 1);
    BOOST_CHECK_EQUAL(y[2], 1.5);
}

BOOST_AUTO_TEST_CASE( cv_sub_scalar )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    
    sub(y, x, 2);

    BOOST_CHECK_EQUAL(y[0], -1);
    BOOST_CHECK_EQUAL(y[1], 0);
    BOOST_CHECK_EQUAL(y[2], 1);
}

BOOST_AUTO_TEST_CASE( cv_add_scalar )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    
    add(y, x, 2);

    BOOST_CHECK_EQUAL(y[0], 3);
    BOOST_CHECK_EQUAL(y[1], 4);
    BOOST_CHECK_EQUAL(y[2], 5);
}

BOOST_AUTO_TEST_CASE( cv_dot )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    float z = dot(x, y);

    BOOST_CHECK_EQUAL(z, 32);
}

BOOST_AUTO_TEST_CASE( cv_add_scaled )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;

    std::vector z(3);
    
    addScaled(z, x, 2, y);

    BOOST_CHECK_EQUAL(z[0], 9);
    BOOST_CHECK_EQUAL(z[1], 12);
    BOOST_CHECK_EQUAL(z[2], 15);
}

BOOST_AUTO_TEST_CASE( cv_square_sum )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    float z = squareSum(x);

    BOOST_CHECK_EQUAL(z, 14);
}

BOOST_AUTO_TEST_CASE( cv_copy )
{
    std::vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    std::vector y(3);
    
    copy(y, x);

    BOOST_CHECK_EQUAL(y[0], 1);
    BOOST_CHECK_EQUAL(y[1], 2);
    BOOST_CHECK_EQUAL(y[2], 3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( current_density_tests )

BOOST_AUTO_TEST_CASE( current_density_dirichlet_boundary_conditions )
{
    MatCOOSymmetric<double> mat(3);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    mat(1, 1) = 4;
    mat(1, 2) = 5;
    mat(2, 2) = 6;

    // 1, 2, 3
    // 2, 4, 5
    // 3, 5, 6

    std::vector b(3);
    b[0] = 14;
    b[1] = 32;
    b[2] = 54;

    CurrentDensitySystem system(mat, b);

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

    std::vector x(3);
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;

    preconditionedSSORConjugateGradientSolver(csr, system.b, x, 1.5, 1e-10, 1000);

    // wolfram alpha: solve {{1, 0, 0}, {0, 1, 0}, {0, 0, 6}} * {x, y, z} = {1, 2, 41}

    // 1, 0, 0
    // 0, 4, 5
    // 0, 5, 6

    // solution: {x, y, z} = {1, 75, -54}

    // check result inside tolerance
    BOOST_CHECK_SMALL(x[0] - 1, 1e-6);
    BOOST_CHECK_SMALL(x[1] - 2, 1e-6);
    BOOST_CHECK_SMALL(x[2] - 41.0/6.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()