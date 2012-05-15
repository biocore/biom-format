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

uint32_t SparseMatFloat::length() {
    return hash.size();
}

items_float SparseMatFloat::keys() {
    items_float all_keys;
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    all_keys.rows = new uint32_t[n_keys];
    all_keys.cols = new uint32_t[n_keys];
    
    for(std::tr1::unordered_map<uint64_t,double>::iterator i = hash.begin(); i != hash.end(); i++) {
        all_keys.rows[count] = DECODE_KEY_ROW(i->first);
        all_keys.cols[count] = DECODE_KEY_COL(i->first);
        count += 1;
    }

    return all_keys;
}

items_float SparseMatFloat::items() {
    items_float all_keys;
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    all_keys.rows = new uint32_t[n_keys];
    all_keys.cols = new uint32_t[n_keys];
    all_keys.values = new double[n_keys];
    
    for(std::tr1::unordered_map<uint64_t,double>::iterator i = hash.begin(); i != hash.end(); i++) {
        all_keys.rows[count] = DECODE_KEY_ROW(i->first);
        all_keys.cols[count] = DECODE_KEY_COL(i->first);
        all_keys.values[count] = i->second;
        count += 1;
    }

    return all_keys;
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

uint32_t SparseMatInt::length() {
    return hash.size();
}

items_int SparseMatInt::keys() {
    items_int all_keys;
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    all_keys.rows = new uint32_t[n_keys];
    all_keys.cols = new uint32_t[n_keys];
    
    for(std::tr1::unordered_map<uint64_t,int32_t>::iterator i = hash.begin(); i != hash.end(); i++) {
        all_keys.rows[count] = DECODE_KEY_ROW(i->first);
        all_keys.cols[count] = DECODE_KEY_COL(i->first);
        count += 1;
    }

    return all_keys;
}

items_int SparseMatInt::items() {
    items_int all_keys;
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    all_keys.rows = new uint32_t[n_keys];
    all_keys.cols = new uint32_t[n_keys];
    all_keys.values = new int32_t[n_keys];
    
    for(std::tr1::unordered_map<uint64_t,int32_t>::iterator i = hash.begin(); i != hash.end(); i++) {
        all_keys.rows[count] = DECODE_KEY_ROW(i->first);
        all_keys.cols[count] = DECODE_KEY_COL(i->first);
        all_keys.values[count] = i->second;
        count += 1;
    }

    return all_keys;
}