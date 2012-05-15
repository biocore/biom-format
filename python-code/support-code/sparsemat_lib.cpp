#include <iostream>
#include "sparsemat_lib.h"

using namespace sparsemat;

// SparseMatFloat
SparseMatFloat::SparseMatFloat() {}

void SparseMatFloat::insert(uint32_t row, uint32_t col, double value) {
    current_key = MAKE_KEY(row,col);
#ifdef _SPARSEMAT_LIB_DEBUG
    std::cout << "Row: " << row << "\tCol: " << col << std::endl;
    std::cout << "Current key: " << current_key << "\tvalue: " << value << std::endl;
#endif
    hash[current_key] = value;
}

double SparseMatFloat::get(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
#ifdef _SPARSEMAT_LIB_DEBUG
    std::cout << "Row: " << row << "\tCol: " << col << std::endl;
    std::cout << "Current key: " << current_key << std::endl;
#endif
    return hash[current_key];
}

void SparseMatFloat::erase(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
    hash.erase(current_key);
}

int SparseMatFloat::contains(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
    std::tr1::unordered_map<uint64_t,double>::const_iterator got = hash.find(current_key);
    
    if(got == hash.end())
        return 0;
    else
        return 1;
}

// SparseMatInt
SparseMatInt::SparseMatInt() {}

void SparseMatInt::insert(uint32_t row, uint32_t col, int32_t value) {
    current_key = MAKE_KEY(row,col);
#ifdef _SPARSEMAT_LIB_DEBUG
    std::cout << "max_size = " << hash.max_size() << std::endl;
    std::cout << "max_bucket_count = " << hash.max_bucket_count() << std::endl;
    std::cout << "max_load_factor = " << hash.max_load_factor() << std::endl;
    std::cout << "Row: " << row << "\tCol: " << col << std::endl;
    std::cout << "Current key: " << current_key << "\tvalue: " << value << std::endl;
#endif
    hash[current_key] = value;
}

int32_t SparseMatInt::get(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
#ifdef _SPARSEMAT_LIB_DEBUG
    std::cout << "Row: " << row << "\tCol: " << col << std::endl;
    std::cout << "Current key: " << current_key << std::endl;
#endif
    return hash[current_key];
}

void SparseMatInt::erase(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
    hash.erase(current_key);
}

int SparseMatInt::contains(uint32_t row, uint32_t col) {
    current_key = MAKE_KEY(row, col);
    std::tr1::unordered_map<uint64_t,int32_t>::const_iterator got = hash.find(current_key);
    
    if(got == hash.end())
        return 0;
    else
        return 1;
}
