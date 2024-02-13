#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <numeric>

#include "math.hpp"

namespace nikfemm {
    void mult(std::vector<double>& result, const MatCSRSymmetric& mat, const std::vector<double>& cv) {
        for (uint32_t i = mat.m; i-- > 0;) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[i] += mat.val[j] * cv[mat.col_ind[j]];
                result[mat.col_ind[j]] += mat.val[j] * cv[i];
            }
        }

        /*
        for (uint32_t i = 0; i < mat.m; i++) {
            result[i] = mat.diag[i] * cv[i];
        }

        for (uint32_t i = 0; i < mat.m; i++) {
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[i] += mat.val[j] * cv[mat.col_ind[j]];
            }
        }

        for (uint32_t i = 0; i < mat.m; i++) {
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[mat.col_ind[j]] += mat.val[j] * cv[i];
            }
        }
        */
    }

    void mult(std::vector<double>& result, const MatCSRLowerTri& mat, const std::vector<double>& cv) {
        for (uint32_t i = mat.m; i-- > 0;) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[mat.col_ind[j]] += mat.val[j] * cv[i];
            }
        }
    }

    void mult(std::vector<double>& result, const MatCSRUpperTri& mat, const std::vector<double>& cv) {
        for (uint32_t i = mat.m; i-- > 0;) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[i] += mat.val[j] * cv[mat.col_ind[j]];
            }
        }
    }

    void mult(std::vector<double>& result, const double d, const std::vector<double>& cv) {
        #ifdef NIK_PARALLEL
        for (size_t i = 0; i < cv.size(); i += 8) {
            result[i] = d * cv[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = d * cv[i];
        }
        #endif
    }

    void mult(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2) {
        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] * cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] * cv2[i];
        }
        #endif
    }

    void add(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2) {
        // assert(cv1.size() == cv2.size());
        // assert(cv1.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] + cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] + cv2[i];
        }
        #endif
    }

    void sub(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2) {
        // assert(cv1.size() == cv2.size());
        // assert(cv1.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] - cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] - cv2[i];
        }
        #endif
    }

    void div(std::vector<double>& result, const std::vector<double>& cv, const double d) {
        // assert(cv.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i] / d;
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i] / d;
        }
        #endif
    }

    void div(std::vector<double>& result, const std::vector<double>& cv1, const std::vector<double>& cv2) {
        // assert(cv1.size() == cv2.size());
        // assert(cv1.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] / cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] / cv2[i];
        }
        #endif
    }
    
    void sub(std::vector<double>& result, const std::vector<double>& cv, const double d) {
        // assert(cv.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i] - d;
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i] - d;
        }
        #endif
    }

    void add(std::vector<double>& result, const std::vector<double>& cv, const double d) {
        // assert(cv.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i] + d;
        }
        #else
        // Original non-parallel code goes here
        // std::transform(cv.val.begin(), cv.val.end(), result.val.begin(), std::bind2nd(std::plus<double>(), d));
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i] + d;
        }
        #endif
    }

    void copy(std::vector<double>& result, const std::vector<double>& cv) {
        // assert(cv.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i];
        }
        #else
        // Original non-parallel code goes here
        std::copy(cv.begin(), cv.end(), result.begin());
        #endif
    }

    double dot(const std::vector<double>& cv1, const std::vector<double>& cv2) {
        // assert(cv1.size() == cv2.size());

        #ifdef NIK_PARALLEL
        double result = 0;
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result += cv1[i] * cv2[i];
        }
        return result;
        #else
        // Original non-parallel code goes here
        return std::inner_product(cv1.begin(), cv1.end(), cv2.begin(), 0.0);
        #endif
    }

    void addScaled(std::vector<double>& result, const std::vector<double>& cv1, const double d, const std::vector<double>& cv2) {
        // assert(cv1.size() == cv2.size());
        // assert(cv1.size() == result.size());

        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] + d * cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.size(); i++) {
            result[i] = cv1[i] + d * cv2[i];
        }
        #endif
    }

    double squareSum(const std::vector<double>& cv) {
        double result = 0;
        #ifdef NIK_PARALLEL
        for (uint32_t i = 0; i < cv.size(); i++) {
            result += cv[i] * cv[i];
        }
        #else
        // Original non-parallel code goes here
        result = std::inner_product(cv.begin(), cv.end(), cv.begin(), 0.0);
        #endif
        return result;
    }

    double norm(const std::vector<double>& cv) {
        return sqrt(squareSum(cv));
    }

    std::vector<double> element_wise_norm(std::vector<Vector>& cv) {
        std::vector<double> result(cv.size());
        for (uint32_t i = 0; i < cv.size(); i++) {
            result[i] = cv[i].norm();
        }
        return result;
    }

    bool isnan(const std::vector<double>& cv) {
        for (uint32_t i = 0; i < cv.size(); i++) {
            if (std::isnan(cv[i])) {
                return true;
            }
        }
        return false;
    }

    bool isnan(const MatCSRSymmetric& mat) {
        for (uint32_t i = 0; i < mat.val.size(); i++) {
            if (std::isnan(mat.val[i])) {
                return true;
            }
        }
        for (uint32_t i = 0; i < mat.diag.size(); i++) {
            if (std::isnan(mat.diag[i])) {
                return true;
            }
        }
        return false;
    }

    bool isnan(const double val) {
        return std::isnan(val);
    }
}