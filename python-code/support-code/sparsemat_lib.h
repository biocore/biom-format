#include <tr1/unordered_map>

#define MAKE_KEY(row,col) (uint64_t(row) << 32 | col)

namespace sparsemat {
    class SparseMatFloat {
        private:
            std::tr1::unordered_map<uint64_t, double> hash;
            uint64_t current_key;
        public:
            SparseMatFloat();
            void insert(uint32_t row, uint32_t col, double value);
            double get(uint32_t row, uint32_t col);
            void erase(uint32_t row, uint32_t col);
            //bool exists(uint32_t row, uint32_t col);
    };

    class SparseMatInt {
        private:
            std::tr1::unordered_map<uint64_t, int> hash;
            uint64_t current_key;
        public:
            SparseMatInt();
            void insert(uint32_t row, uint32_t col, int value);
            int get(uint32_t row, uint32_t col);
            void erase(uint32_t row, uint32_t col);
            //bool exists(uint32_t row, uint32_t col);
    };
}
