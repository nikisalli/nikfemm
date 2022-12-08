#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <numeric>

#include <omp.h>

#include "simple_vector.hpp"

namespace nikfemm {
    CV::CV(uint32_t size) {
        val = std::vector<double>(size);
    }

    CV::CV() {
        val = std::vector<double>();
    }

    CV::~CV() {

    }

    void CV::print() const {
        printf("[");
        for (uint32_t i = 0; i < val.size(); i++) {
            printf("%.17g ", val[i]);
        }
        printf("]");
        printf("\n");
    }

    void CV::write_to_file(const char *filename) {
        FILE *f = fopen(filename, "w");
        for (uint32_t i = 0; i < val.size(); i++) {
            fprintf(f, "%.17g ", val[i]);
        }
        fclose(f);
    }

    void CV::mult(CV& result, const MatCSRSymmetric& mat, const CV& cv) {
        for (uint32_t i = mat.m; i-- > 0;) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[i] += mat.val[j] * cv[mat.col_ind[j]];
                result[mat.col_ind[j]] += mat.val[j] * cv[i];
            }
        }
    }

    void CV::mult(CV& result, const MatCSRLowerTri& mat, const CV& cv) {
        for (uint32_t i = mat.m; i-- > 0;) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[mat.col_ind[j]] += mat.val[j] * cv[i];
            }
        }
    }

    void CV::mult(CV& result, const MatCSRUpperTri& mat, const CV& cv) {
        for (uint32_t i = mat.m; i-- > 0;) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[i] += mat.val[j] * cv[mat.col_ind[j]];
            }
        }
    }

    void CV::mult(CV& result, const double d, const CV& cv) {
        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = d * cv[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = d * cv[i];
        }
        #endif
    }


    void CV::mult(CV& result, const CV& cv1, const CV& cv2) {
        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] * cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] * cv2[i];
        }
        #endif
    }

    void CV::add(CV& result, const CV& cv1, const CV& cv2) {
        // assert(cv1.val.size() == cv2.val.size());
        // assert(cv1.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] + cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] + cv2[i];
        }
        #endif
    }

    void CV::sub(CV& result, const CV& cv1, const CV& cv2) {
        // assert(cv1.val.size() == cv2.val.size());
        // assert(cv1.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] - cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] - cv2[i];
        }
        #endif
    }

    void CV::div(CV& result, const CV& cv, const double d) {
        // assert(cv.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i] / d;
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i] / d;
        }
        #endif
    }

    void CV::div(CV& result, const CV& cv1, const CV& cv2) {
        // assert(cv1.val.size() == cv2.val.size());
        // assert(cv1.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] / cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] / cv2[i];
        }
        #endif
    }
    
    void CV::sub(CV& result, const CV& cv, const double d) {
        // assert(cv.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i] - d;
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i] - d;
        }
        #endif
    }

    void CV::add(CV& result, const CV& cv, const double d) {
        // assert(cv.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i] + d;
        }
        #else
        // Original non-parallel code goes here
        // std::transform(cv.val.begin(), cv.val.end(), result.val.begin(), std::bind2nd(std::plus<double>(), d));
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i] + d;
        }
        #endif
    }

    void CV::copy(CV& result, const CV& cv) {
        // assert(cv.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result[i] = cv[i];
        }
        #else
        // Original non-parallel code goes here
        std::copy(cv.val.begin(), cv.val.end(), result.val.begin());
        #endif
    }

    double CV::dot(const CV& cv1, const CV& cv2) {
        // assert(cv1.val.size() == cv2.val.size());

        #ifdef NIK_PARALLEL
        double result = 0;
        #pragma omp parallel for reduction(+:result)
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result += cv1[i] * cv2[i];
        }
        return result;
        #else
        // Original non-parallel code goes here
        return std::inner_product(cv1.val.begin(), cv1.val.end(), cv2.val.begin(), 0.0);
        #endif
    }

    void CV::addScaled(CV& result, const CV& cv1, const double d, const CV& cv2) {
        // assert(cv1.val.size() == cv2.val.size());
        // assert(cv1.val.size() == result.val.size());

        #ifdef NIK_PARALLEL
        #pragma omp parallel for
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] + d * cv2[i];
        }
        #else
        // Original non-parallel code goes here
        for (uint32_t i = 0; i < cv1.val.size(); i++) {
            result[i] = cv1[i] + d * cv2[i];
        }
        #endif
    }

    double CV::squareSum(const CV& cv) {
        double result = 0;
        #ifdef NIK_PARALLEL
        #pragma omp parallel for reduction(+:result)
        for (uint32_t i = 0; i < cv.val.size(); i++) {
            result += cv[i] * cv[i];
        }
        #else
        // Original non-parallel code goes here
        result =  std::inner_product(cv.val.begin(), cv.val.end(), cv.val.begin(), 0.0);
        #endif
        return result;
    }

    double CV::norm(const CV& cv) {
        return sqrt(CV::squareSum(cv));
    }

    void CV::add_elem(uint32_t _m, double d) {
        if (_m > val.size()) {
            val.resize(_m);
        }
        this->val[_m] += d;
    }

    void CV::set_elem(uint32_t _m, double d) {
        if (_m > val.size()) {
            val.resize(_m);
        }
        this->val[_m] = d;
    }
}