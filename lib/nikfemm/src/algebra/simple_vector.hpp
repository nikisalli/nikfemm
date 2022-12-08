#ifndef NIK_COLUMN_VECTOR_HPP
#define NIK_COLUMN_VECTOR_HPP

#include <cstdint>
#include <vector>

#include "csr.hpp"

namespace nikfemm {
    struct MatCSRSymmetric;
    struct MatCSRLowerTri;
    struct MatCSRUpperTri;

    struct CV {
        protected:

        public:
        std::vector<double> val;
        CV(uint32_t size);
        CV();
        ~CV();

        void print() const;
        void write_to_file(const char *filename);
        
        inline double& operator[](uint32_t i) {
            return val[i];
        }

        inline double operator[](uint32_t i) const {
            return val[i];
        }

        static void mult(CV& result, const MatCSRSymmetric& mat, const CV& cv);
        static void mult(CV& result, const MatCSRLowerTri& mat, const CV& cv);
        static void mult(CV& result, const MatCSRUpperTri& mat, const CV& cv);
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
        static double norm(const CV& cv);

        void add_elem(uint32_t _m, double d);  // row, value
        void set_elem(uint32_t _m, double d);  // row, value
    };
}

#endif