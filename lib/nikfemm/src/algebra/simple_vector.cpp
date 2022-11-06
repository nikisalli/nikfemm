#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "simple_vector.hpp"

namespace nikfemm {
    CV::CV(uint32_t size) {
        val = new double[size]();
        m = size;
    }

    CV::~CV() {

    }

    void CV::print() {
        printf("[");
        for (uint32_t i = 0; i < m; i++) {
            printf("%.17g ", val[i]);
        }
        printf("]");
        printf("\n");
    }

    void CV::write_to_file(const char *filename) {
        FILE *f = fopen(filename, "w");
        for (uint32_t i = 0; i < m; i++) {
            fprintf(f, "%.17g ", val[i]);
        }
        fclose(f);
    }

    double& CV::operator[](uint32_t i) {
        return val[i];
    }

    double CV::operator[](uint32_t i) const {
        return val[i];
    }

    void CV::mult(CV& result, const MatCSR& mat, const CV& cv) {
        for (uint32_t i = 0; i < mat.m; i++) {
            result[i] = 0;
            for (uint32_t j = mat.IA[i]; j < mat.IA[i + 1]; j++) {
                result[i] += mat.A[j] * cv[mat.JA[j]];
            }
        }
    }
    
    void CV::mult(CV& result, const MatSSS& mat, const CV& cv) {
        // SSS is symmetric, so we can use the lower triangular part only
        for (uint32_t i = 0; i < mat.m; i++) {
            result[i] = mat.diag[i] * cv[i];
            for (uint32_t j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; j++) {
                result[i] += mat.val[j] * cv[mat.col_ind[j]];
                result[mat.col_ind[j]] += mat.val[j] * cv[i];
            }
        }
    }

    void CV::mult(CV& result, const double d, const CV& cv) {
        for (uint32_t i = 0; i < cv.m; i++) {
            result[i] = d * cv[i];
        }
    }

    void CV::mult(CV& result, const CV& cv1, const CV& cv2) {
        for (uint32_t i = 0; i < cv1.m; i++) {
            result[i] = cv1[i] * cv2[i];
        }
    }

    void CV::add(CV& result, const CV& cv1, const CV& cv2) {
        assert(cv1.m == cv2.m);
        assert(cv1.m == result.m);
        for (uint32_t i = 0; i < cv1.m; i++) {
            result[i] = cv1[i] + cv2[i];
        }
    }

    void CV::sub(CV& result, const CV& cv1, const CV& cv2) {
        assert(cv1.m == cv2.m);
        assert(cv1.m == result.m);
        for (uint32_t i = 0; i < cv1.m; i++) {
            result[i] = cv1[i] - cv2[i];
        }
    }

    void CV::div(CV& result, const CV& cv, const double d) {
        assert(cv.m == result.m);
        for (uint32_t i = 0; i < cv.m; i++) {
            result[i] = cv[i] / d;
        }
    }

    void CV::div(CV& result, const CV& cv1, const CV& cv2) {
        assert(cv1.m == cv2.m);
        assert(cv1.m == result.m);
        for (uint32_t i = 0; i < cv1.m; i++) {
            result[i] = cv1[i] / cv2[i];
        }
    }
    
    void CV::sub(CV& result, const CV& cv, const double d) {
        assert(cv.m == result.m);
        for (uint32_t i = 0; i < cv.m; i++) {
            result[i] = cv[i] - d;
        }
    }

    void CV::add(CV& result, const CV& cv, const double d) {
        assert(cv.m == result.m);
        for (uint32_t i = 0; i < cv.m; i++) {
            result[i] = cv[i] + d;
        }
    }

    void CV::copy(CV& result, const CV& cv) {
        assert(cv.m == result.m);
        for (uint32_t i = 0; i < cv.m; i++) {
            result[i] = cv[i];
        }
    }

    double CV::dot(const CV& cv1, const CV& cv2) {
        assert(cv1.m == cv2.m);
        double result = 0;
        for (uint32_t i = 0; i < cv1.m; i++) {
            result += cv1[i] * cv2[i];
        }
        return result;
    }

    void CV::addScaled(CV& result, const CV& cv1, const double d, const CV& cv2) {
        assert(cv1.m == cv2.m);
        assert(cv1.m == result.m);
        for (uint32_t i = 0; i < cv1.m; i++) {
            result[i] = cv1[i] + d * cv2[i];
        }
    }

    double CV::squareSum(const CV& cv) {
        double result = 0;
        for (uint32_t i = 0; i < cv.m; i++) {
            result += cv[i] * cv[i];
        }
        return result;
    }

    void CV::add_elem(uint32_t _m, double d) {
        if (_m > m) m = _m;
        this->val[_m] += d;
    }

    void CV::set_elem(uint32_t _m, double d) {
        if (_m > m) m = _m;
        this->val[_m] = d;
    }
}