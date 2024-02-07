#ifndef NIKFEMM_ALGEBRA_MATH_HPP
#define NIKFEMM_ALGEBRA_MATH_HPP

#include <vector>

#include "csr.hpp"

namespace nikfemm {
    void mult(std::vector<double>& result, const MatCSRSymmetric& mat, const std::vector<double>& cv);
    void mult(std::vector<double>& result, const MatCSRLowerTri& mat, const std::vector<double>& cv);
    void mult(std::vector<double>& result, const MatCSRUpperTri& mat, const std::vector<double>& cv);
    void mult(std::vector<double>& result, const double d, const std::vector<double>& cv);
    void mult(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2);
    void add(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2);
    void sub(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2);
    void div(std::vector<double>& result, const std::vector<double>& cv, const double d);
    void div(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2);
    void sub(std::vector<double>& result, const std::vector<double>& cv, const double d);
    void add(std::vector<double>& result, const std::vector<double>& cv, const double d);
    void addScaled(std::vector<double>& result, const std::vector<double>& cv1, const double d, const std::vector<double>& cv2);
    void copy(std::vector<double>& result, const std::vector<double>& cv);

    double squareSum(const std::vector<double>& cv);
    double dot(const std::vector<double>& cv1, const std::vector<double>& cv2);
    double norm(const std::vector<double>& cv);
}

#endif // NIKFEMM_ALGEBRA_MATH_HPP