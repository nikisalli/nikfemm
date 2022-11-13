#include <nikfemm.hpp>

#define BOOST_TEST_MODULE VsidCommonTest
#include <boost/test/unit_test.hpp>

using namespace nikfemm;

BOOST_AUTO_TEST_SUITE( csr_tests )

BOOST_AUTO_TEST_CASE( csr_index )
{
    MatCOO coo(3);
    coo.set_elem(0, 0, 1);
    coo.set_elem(0, 1, 2);
    coo.set_elem(1, 0, 2);
    coo.set_elem(1, 2, 4);
    coo.set_elem(2, 1, 4);
    coo.set_elem(1, 1, 5);
    coo.set_elem(2, 2, 6);

    MatCSR csr(coo);

    coo.print();
    csr.print();
    csr.printCSR();

    BOOST_CHECK_EQUAL(csr(0, 0), 1);
    BOOST_CHECK_EQUAL(csr(0, 1), 2);
    BOOST_CHECK_EQUAL(csr(1, 0), 2);
    BOOST_CHECK_EQUAL(csr(1, 2), 4);
    BOOST_CHECK_EQUAL(csr(2, 1), 4);
    BOOST_CHECK_EQUAL(csr(1, 1), 5);
    BOOST_CHECK_EQUAL(csr(2, 2), 6);
}

BOOST_AUTO_TEST_CASE( csr_conjugate_gradient_solve )
{
    MatCOO coo(2);
    // [1, 2, 0; 0, 3, 4; 5, 0, 6]
    coo.set_elem(0, 0, 3);
    coo.set_elem(0, 1, 1);
    coo.set_elem(1, 0, 1);
    coo.set_elem(1, 1, 2);

    CV b(2);
    b[0] = -1;
    b[1] = 2;

    CV x0(2);
    x0[0] = 0;
    x0[1] = 0;

    MatCSR csr(coo);

    csr.conjugateGradientSolver(b, x0, 1e-6, 100);

    x0.print();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( sss_tests )

BOOST_AUTO_TEST_CASE( sss_index ) {
    MatCOO coo(8);
    // [1, 2, 0; 0, 3, 4; 5, 0, 6]
    coo.set_elem(0, 0, 1);
    coo.set_elem(1, 1, 2);
    coo.set_elem(2, 2, 3);
    coo.set_elem(3, 3, 4);
    coo.set_elem(4, 4, 5);
    coo.set_elem(5, 5, 6);
    coo.set_elem(6, 6, 7);
    coo.set_elem(7, 7, 8);
    coo.set_elem(2, 0, 9);
    coo.set_elem(6, 0, 10);
    coo.set_elem(5, 2, 11);
    coo.set_elem(7, 4, 12);

    MatSSS sss(coo);
    coo.print();

    sss.printSSS();
    sss.print();

    BOOST_CHECK_EQUAL(sss(0, 0), 1);
    BOOST_CHECK_EQUAL(sss(0, 1), 0);
    BOOST_CHECK_EQUAL(sss(0, 2), 9);
    BOOST_CHECK_EQUAL(sss(0, 3), 0);
    BOOST_CHECK_EQUAL(sss(0, 4), 0);
    BOOST_CHECK_EQUAL(sss(0, 5), 0);
    BOOST_CHECK_EQUAL(sss(0, 6), 10);
    BOOST_CHECK_EQUAL(sss(0, 7), 0);
    BOOST_CHECK_EQUAL(sss(1, 0), 0);
    BOOST_CHECK_EQUAL(sss(1, 1), 2);
    BOOST_CHECK_EQUAL(sss(1, 2), 0);
    BOOST_CHECK_EQUAL(sss(1, 3), 0);
    BOOST_CHECK_EQUAL(sss(1, 4), 0);
    BOOST_CHECK_EQUAL(sss(1, 5), 0);
    BOOST_CHECK_EQUAL(sss(1, 6), 0);
    BOOST_CHECK_EQUAL(sss(1, 7), 0);
    BOOST_CHECK_EQUAL(sss(2, 0), 9);
    BOOST_CHECK_EQUAL(sss(2, 1), 0);
    BOOST_CHECK_EQUAL(sss(2, 2), 3);
    BOOST_CHECK_EQUAL(sss(2, 3), 0);
    BOOST_CHECK_EQUAL(sss(2, 4), 0);
    BOOST_CHECK_EQUAL(sss(2, 5), 11);
    BOOST_CHECK_EQUAL(sss(2, 6), 0);
    BOOST_CHECK_EQUAL(sss(2, 7), 0);
    BOOST_CHECK_EQUAL(sss(3, 0), 0);
    BOOST_CHECK_EQUAL(sss(3, 1), 0);
    BOOST_CHECK_EQUAL(sss(3, 2), 0);
    BOOST_CHECK_EQUAL(sss(3, 3), 4);
    BOOST_CHECK_EQUAL(sss(3, 4), 0);
    BOOST_CHECK_EQUAL(sss(3, 5), 0);
    BOOST_CHECK_EQUAL(sss(3, 6), 0);
    BOOST_CHECK_EQUAL(sss(3, 7), 0);
    BOOST_CHECK_EQUAL(sss(4, 0), 0);
    BOOST_CHECK_EQUAL(sss(4, 1), 0);
    BOOST_CHECK_EQUAL(sss(4, 2), 0);
    BOOST_CHECK_EQUAL(sss(4, 3), 0);
    BOOST_CHECK_EQUAL(sss(4, 4), 5);
    BOOST_CHECK_EQUAL(sss(4, 5), 0);
    BOOST_CHECK_EQUAL(sss(4, 6), 0);
    BOOST_CHECK_EQUAL(sss(4, 7), 12);
    BOOST_CHECK_EQUAL(sss(5, 0), 0);
    BOOST_CHECK_EQUAL(sss(5, 1), 0);
    BOOST_CHECK_EQUAL(sss(5, 2), 11);
    BOOST_CHECK_EQUAL(sss(5, 3), 0);
    BOOST_CHECK_EQUAL(sss(5, 4), 0);
    BOOST_CHECK_EQUAL(sss(5, 5), 6);
    BOOST_CHECK_EQUAL(sss(5, 6), 0);
    BOOST_CHECK_EQUAL(sss(5, 7), 0);
    BOOST_CHECK_EQUAL(sss(6, 0), 10);
    BOOST_CHECK_EQUAL(sss(6, 1), 0);
    BOOST_CHECK_EQUAL(sss(6, 2), 0);
    BOOST_CHECK_EQUAL(sss(6, 3), 0);
    BOOST_CHECK_EQUAL(sss(6, 4), 0);
    BOOST_CHECK_EQUAL(sss(6, 5), 0);
    BOOST_CHECK_EQUAL(sss(6, 6), 7);
    BOOST_CHECK_EQUAL(sss(6, 7), 0);
    BOOST_CHECK_EQUAL(sss(7, 0), 0);
    BOOST_CHECK_EQUAL(sss(7, 1), 0);
    BOOST_CHECK_EQUAL(sss(7, 2), 0);
    BOOST_CHECK_EQUAL(sss(7, 3), 0);
    BOOST_CHECK_EQUAL(sss(7, 4), 12);
    BOOST_CHECK_EQUAL(sss(7, 5), 0);
    BOOST_CHECK_EQUAL(sss(7, 6), 0);
    BOOST_CHECK_EQUAL(sss(7, 7), 8);
}

BOOST_AUTO_TEST_CASE( sss_cv_mult ) {
    MatCOO coo(8);
    // [1, 2, 0; 0, 3, 4; 5, 0, 6]
    coo.set_elem(0, 0, 1);
    coo.set_elem(1, 1, 2);
    coo.set_elem(2, 2, 3);
    coo.set_elem(3, 3, 4);
    coo.set_elem(4, 4, 5);
    coo.set_elem(5, 5, 6);
    coo.set_elem(6, 6, 7);
    coo.set_elem(7, 7, 8);

    coo.set_elem(2, 0, 9);
    coo.set_elem(6, 0, 10);
    coo.set_elem(5, 2, 11);
    coo.set_elem(7, 4, 12);

    coo.set_elem(0, 2, 9);
    coo.set_elem(0, 6, 10);
    coo.set_elem(2, 5, 11);
    coo.set_elem(4, 7, 12);

    MatSSS sss(coo);

    sss.print();

    CV cv(8);
    cv.set_elem(0, 1);
    cv.set_elem(1, 2);
    cv.set_elem(2, 3);
    cv.set_elem(3, 4);
    cv.set_elem(4, 5);
    cv.set_elem(5, 6);
    cv.set_elem(6, 7);
    cv.set_elem(7, 8);

    CV cv_res(8);

    CV::mult(cv_res, sss, cv);
    
    cv_res.print();
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
    MatCOO coo(3);
    coo.set_elem(0, 0, 1);
    coo.set_elem(0, 1, 2);
    coo.set_elem(1, 0, 2);
    coo.set_elem(1, 2, 4);
    coo.set_elem(2, 1, 4);
    coo.set_elem(1, 1, 5);
    coo.set_elem(2, 2, 6);

    MatCSR csr(coo);

    CV x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    CV y(3);
    
    CV::mult(y, csr, x);

    BOOST_CHECK_EQUAL(y[0], 5);
    BOOST_CHECK_EQUAL(y[1], 24);
    BOOST_CHECK_EQUAL(y[2], 26);
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