#include <iostream>
#include "sparsemat_lib.h"
#include <stdio.h>
#include <vector>
using namespace sparsemat;

// SparseMatFloat
SparseMatFloat::SparseMatFloat() {
    tmp_items.rows = NULL;
    tmp_items.cols = NULL;
    tmp_items.values = NULL;
}

SparseMatFloat::~SparseMatFloat() {
    cleanTmpItems();
}

void SparseMatFloat::insert(uint32_t row, uint32_t col, double value) {
    if(value == 0.0)
        return;
        
    current_key = MAKE_KEY(row,col);
    hash[current_key] = value;
}

double SparseMatFloat::get(uint32_t row, uint32_t col) {
    // have to short circuit, the [] on return adds the item. WTF
    if(!contains(row,col))
        return 0.0;
        
    current_key = MAKE_KEY(row, col);
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
    uint32_t n_keys = length();
    uint32_t count = 0;
        
    cleanTmpItems();
    
    tmp_items.rows = new uint32_t[n_keys];
    tmp_items.cols = new uint32_t[n_keys];
    tmp_items.values = NULL;
    
    for(std::tr1::unordered_map<uint64_t,double>::iterator i = hash.begin(); i != hash.end(); i++) {
        tmp_items.rows[count] = DECODE_KEY_ROW(i->first);
        tmp_items.cols[count] = DECODE_KEY_COL(i->first);
    
        count += 1;
    }

    return tmp_items;
}

items_float SparseMatFloat::items() {
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    cleanTmpItems();
    
    tmp_items.rows = new uint32_t[n_keys];
    tmp_items.cols = new uint32_t[n_keys];
    tmp_items.values = new double[n_keys];
    
    for(std::tr1::unordered_map<uint64_t,double>::iterator i = hash.begin(); i != hash.end(); i++) {
        tmp_items.rows[count] = DECODE_KEY_ROW(i->first);
        tmp_items.cols[count] = DECODE_KEY_COL(i->first);
        tmp_items.values[count] = i->second;
        count += 1;
    }

    return tmp_items;
}

void SparseMatFloat::cleanTmpItems() {
    if(tmp_items.rows != NULL)
        delete tmp_items.rows;
    
    if(tmp_items.cols != NULL)
        delete tmp_items.cols;
    
    if(tmp_items.values != NULL)
        delete tmp_items.values;
}

// SparseMatInt
SparseMatInt::SparseMatInt() {
    tmp_items.rows = NULL;
    tmp_items.cols = NULL;
    tmp_items.values = NULL;
}

SparseMatInt::~SparseMatInt() {
    cleanTmpItems();
}

void SparseMatInt::insert(uint32_t row, uint32_t col, int32_t value) {
    if(value == 0)
        return;
        
    current_key = MAKE_KEY(row,col);
    hash[current_key] = value;
}

int32_t SparseMatInt::get(uint32_t row, uint32_t col) {
    // have to short circuit, the [] on return adds the item. WTF
    if(!contains(row,col))
        return 0;

    current_key = MAKE_KEY(row, col);
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
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    cleanTmpItems();
        
    tmp_items.rows = new uint32_t[n_keys];
    tmp_items.cols = new uint32_t[n_keys];
    tmp_items.values = NULL;
    
    for(std::tr1::unordered_map<uint64_t,int32_t>::iterator i = hash.begin(); i != hash.end(); i++) {
        tmp_items.rows[count] = DECODE_KEY_ROW(i->first);
        tmp_items.cols[count] = DECODE_KEY_COL(i->first);
        count += 1;
    }

    // ARG, return by value...
    return tmp_items;
}

items_int SparseMatInt::items() {
    uint32_t n_keys = length();
    uint32_t count = 0;
    
    cleanTmpItems();

    tmp_items.rows = new uint32_t[n_keys];
    tmp_items.cols = new uint32_t[n_keys];
    tmp_items.values = new int32_t[n_keys];
    
    for(std::tr1::unordered_map<uint64_t,int32_t>::iterator i = hash.begin(); i != hash.end(); i++) {
        tmp_items.rows[count] = DECODE_KEY_ROW(i->first);
        tmp_items.cols[count] = DECODE_KEY_COL(i->first);
        tmp_items.values[count] = i->second;
        count += 1;
    }

    return tmp_items;
}

void SparseMatInt::cleanTmpItems() {   
    if(tmp_items.rows != NULL)
        delete tmp_items.rows;
    
    if(tmp_items.cols != NULL)
        delete tmp_items.cols;
    
    if(tmp_items.values != NULL)
        delete tmp_items.values;
}
