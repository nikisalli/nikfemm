#ifndef NIK_COLUMN_VECTOR_HPP
#define NIK_COLUMN_VECTOR_HPP

#include <cstdint>
#include <vector>

#include "csr.hpp"
#include "sss.hpp"

namespace nikfemm {
    struct MatCSR;
    struct MatSSS;
    struct CV {
        double* val;
        uint32_t m;  // columns

        CV(uint32_t size);
        ~CV();

        void print();
        void write_to_file(const char *filename);
        
        double& operator[](uint32_t i);
        double operator[](uint32_t i) const;

        static void mult(CV& result, const MatCSR& mat, const CV& cv);
        static void mult(CV& result, const MatSSS& mat, const CV& cv);
        static void mult(CV& result, const double d, const CV& cv);
        static void mult(CV& result, const CV& cv1, const CV& cv2);
        static void add(CV& result, const CV& cv1, const CV& cv2);
        static void sub(CV& result, const CV& cv1, const CV& cv2);
        static void div(CV& result, const CV& cv, const double d);
        static void div(CV& result, const CV& cv1, const CV& cv2);
        static void sub(CV& result, const CV& cv, const double d);
        static void add(CV& result, const CV& cv, const double d);
        static void addScaled(CV& result, const CV& cv1, const double d, const CV& cv2);
        static void copy(CV& result, const CV& cv);

        static double squareSum(const CV& cv);
        static double dot(const CV& cv1, const CV& cv2);

        void add_elem(uint32_t _m, double d);  // row, value
        void set_elem(uint32_t _m, double d);  // row, value
    };
}

#endif