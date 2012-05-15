#include <iostream>
#include "sparsemat_lib.h"

using namespace sparsemat;

// SparseMatFloat
SparseMatFloat::SparseMatFloat() {}

void SparseMatFloat::insert(uint32_t row, uint32_t col, double value) {
        current_key = MAKE_KEY(row,col);
        hash[current_key] = value;
}

double SparseMatFloat::get(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
    return hash[current_key];
}

// SparseMatInt
SparseMatInt::SparseMatInt() {}

void SparseMatInt::insert(uint32_t row, uint32_t col, int value) {
        current_key = MAKE_KEY(row,col);
        hash[current_key] = value;
}

int SparseMatInt::get(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
    return hash[current_key];
}

