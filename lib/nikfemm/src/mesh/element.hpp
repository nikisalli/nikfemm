#ifndef NIK_ELEMENT_HPP
#define NIK_ELEMENT_HPP

#define EMPTY_ELEMENT nullptr

namespace nikfemm {
    enum ElementType {
        ELEMENT_TYPE_TRIANGLE,
    };

    class Element {
        public:
            double mu_r;

            ElementType type;

            virtual ~Element() = 0;

            virtual double getArea();

            virtual bool operator==(const Element& e) const;
            virtual bool operator!=(const Element& e) const;
    };
}

#endif