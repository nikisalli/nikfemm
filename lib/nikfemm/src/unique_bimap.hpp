#ifndef NIKBIMAP_H
#define NIKBIMAP_H

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <stdexcept>

namespace nikfemm {
    template <typename T1, typename T2>
    class unique_bimap {
        /* bimap with unique keys and values */
        private:
            std::map<T1, T2> map1;
            std::map<T2, T1> map2;
        
        public:
            unique_bimap() {
                map1 = std::map<T1, T2>();
                map2 = std::map<T2, T1>();
            }
            
            ~unique_bimap() {
                map1.clear();
                map2.clear();
            }
            
            void insert(T1 key1, T2 key2) {
                // check if key1 already exists
                if (map1.find(key1) != map1.end()) {
                    fprintf(stderr, "bimap::insert: key1 already exists\n");
                    throw std::invalid_argument("bimap::insert: key1 already exists");
                }
                // check if key2 already exists
                if (map2.find(key2) != map2.end()) {
                    fprintf(stderr, "bimap::insert: key2 already exists\n");
                    throw std::invalid_argument("bimap::insert: key2 already exists");
                }
                map1.insert(std::make_pair(key1, key2));
                map2.insert(std::make_pair(key2, key1));
            }
            
            T2 getLeft(T1 key1) {
                return map1[key1];
            }
            
            T1 getRight(T2 key2) {
                return map2[key2];
            }
            
            void removeLeft(T1 key1) {
                T2 key2 = map1[key1];
                map1.erase(key1);
                map2.erase(key2);
            }
            
            void removeRight(T2 key2) {
                T1 key1 = map2[key2];
                map2.erase(key2);
                map1.erase(key1);
            }
            
            bool containsLeft(T1 key1) {
                return map1.find(key1) != map1.end();
            }
            
            bool containsRight(T2 key2) {
                return map2.find(key2) != map2.end();
            }
            
            void clear() {
                map1.clear();
                map2.clear();
            }
    };
}

#endif // NIKBIMAP_H