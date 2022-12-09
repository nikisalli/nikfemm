#include "coo.hpp"

namespace nikfemm {
    BaseCOO::BaseCOO() {
        elems = std::vector<ElemCOO>();
    }

    BaseCOO::BaseCOO(BuildMatCOO<double>& coo) {
        // assume sorted
        // copy from unordered map to vector
        elems = std::vector<ElemCOO>();
        elems.reserve(coo.elems.size());
        for (auto it = coo.elems.begin(); it != coo.elems.end(); it++) {
            uint64_t key = it->first;
            double value = it->second;
            elems.push_back(ElemCOO(key >> 32, key & 0xFFFFFFFF, value));
        }
        m = coo.m;

        // sort elems by row, col
        std::sort(elems.begin(), elems.end(), [](const auto& a, const auto& b) {
            if (a.row == b.row) {
                return a.row < b.row;
            } else {
                return a.col < b.col;
            }
        });
    }

    BaseCOO::BaseCOO(const BaseCOO& coo) {
        elems = coo.elems;
        m = coo.m;
    }

    BaseCOO::~BaseCOO() {
    }

    void BaseCOO::printCOO() {
        for (uint32_t i = 0; i < elems.size(); i++) {
            printf("(%u, %u, %.1f) ", elems[i].row, elems[i].col, elems[i].val);
        }
    }

    void BaseCOO::print() {
        // assume sorted
        uint32_t i = 0;
        for (uint32_t row = 0; row < m; row++) {
            for (uint32_t col = 0; col < m; col++) {
                if (i < elems.size() && elems[i].row == row && elems[i].col == col) {
                    printf("%.1f ", elems[i].val);
                    i++;
                } else {
                    printf("0 ");
                }
            }
            printf("\n");
        }
    }

    void BaseCOO::write_to_file(const char *filename) {
        FILE *f = fopen(filename, "w");
        // assume sorted
        uint32_t i = 0;
        for (uint32_t row = 0; row < m; row++) {
            for (uint32_t col = 0; col < m; col++) {
                if (i < elems.size() && elems[i].row == row && elems[i].col == col) {
                    fprintf(f, "%.1f ", elems[i].val);
                    i++;
                } else {
                    fprintf(f, "0 ");
                }
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
}