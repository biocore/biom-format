#include <tr1/unordered_map>

#define MAKE_KEY(row,col) (uint64_t(row) << 32 | col)
#define DECODE_KEY_ROW(key) uint32_t(key >> 32)
#define DECODE_KEY_COL(key) uint32_t(key & 0x00000000ffffffff)

namespace sparsemat {
    typedef std::tr1::unordered_map<uint64_t, double> float_map_t;
    typedef std::tr1::unordered_map<uint64_t, int32_t> int_map_t;
    
    struct items_int {
        uint32_t *rows;
        uint32_t *cols;
        int32_t *values;
    };

    struct items_float {
        uint32_t *rows;
        uint32_t *cols;
        double *values;
    };

    class SparseMatFloat {
        private:
            float_map_t hash;
            uint64_t current_key;
        public:
            SparseMatFloat();
            void insert(uint32_t row, uint32_t col, double value);
            double get(uint32_t row, uint32_t col);
            void erase(uint32_t row, uint32_t col);
            int contains(uint32_t row, uint32_t col);
            uint32_t length();
            items_float keys();
    };

    class SparseMatInt {
        private:
            int_map_t hash;
            uint64_t current_key;
        public:
            SparseMatInt();
            void insert(uint32_t row, uint32_t col, int32_t value);
            int32_t get(uint32_t row, uint32_t col);
            void erase(uint32_t row, uint32_t col);
            int contains(uint32_t row, uint32_t col);
            uint32_t length();
            items_int keys();
    };
}
