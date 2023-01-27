#include <nikfemm.hpp>

#define BOOST_TEST_MODULE VsidCommonTest
#include <boost/test/unit_test.hpp>

using namespace nikfemm;

BOOST_AUTO_TEST_SUITE( csr_tests )

BOOST_AUTO_TEST_CASE( csr_index )
{
    
}

BOOST_AUTO_TEST_CASE( csr_cv_mult )
{
   
}

BOOST_AUTO_TEST_CASE( csr_conjugate_gradient_solve )
{
    
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