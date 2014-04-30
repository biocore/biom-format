# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division

from itertools import compress

import numpy as np
cimport numpy as cnp


cdef _zero_rows_csr(arr, cnp.ndarray[cnp.int8_t, ndim=1] booleans):
    cdef Py_ssize_t m, n, row, j
    cdef cnp.int32_t start, end
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    m, n = arr.shape
    for row in range(m):
        if booleans[row]:
            continue
        start, end = indptr[row], indptr[row+1]
        for j in range(start, end):
            data[j] = 0
    arr.eliminate_zeros()
    
cdef _zero_columns_csc(arr, cnp.ndarray[cnp.int8_t, ndim=1] booleans):
    cdef Py_ssize_t m, n, col, j
    cdef cnp.int32_t start, end
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    m, n = arr.shape
    for col in range(n):
        if booleans[col]:
            continue
        start, end = indptr[col], indptr[col+1]
        for j in range(start, end):
            data[j] = 0
    arr.eliminate_zeros()

cdef _zero_columns_csr(arr, cnp.ndarray[cnp.int8_t, ndim=1] booleans):
    cdef Py_ssize_t i, col
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indices = arr.indices
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    for i in range(arr.indices.size):
        col = indices[i]
        if booleans[col]:
            continue
        data[i] = 0
    arr.eliminate_zeros()
    
cdef _zero_rows_csc(arr, cnp.ndarray[cnp.int8_t, ndim=1] booleans):
    cdef Py_ssize_t i, row
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indices = arr.indices
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    for i in range(arr.indices.size):
        row = indices[i]
        if booleans[row]:
            continue
        data[i] = 0
    arr.eliminate_zeros()

cdef cnp.ndarray[cnp.int8_t, ndim=1] _make_filter_array(ids, metadata, func, bint invert):
    cdef cnp.ndarray[cnp.int8_t, ndim=1] bools = np.empty(len(ids), dtype=np.int8)
    if func
    for i in range(len(ids)):
        bools[i] = func(ids[i], metadata[i]) ^ invert
    return bools

cdef _remove_rows_csr(arr, cnp.ndarray[cnp.int8_t, ndim=1] booleans):
    cdef Py_ssize_t m, n, row, j, offset, offset_rows, nnz
    cdef cnp.int32_t start, end
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indices = arr.indices
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    m, n = arr.shape
    offset_rows = 0
    offset = 0
    nnz = 0
    for row in range(m):
        start, end = indptr[row], indptr[row+1]
        if booleans[row]:
            indptr[row-offset_rows] = nnz
            nnz += end - start
            indptr[row-offset_rows + 1] = nnz
            for j in range(start, end):
                data[j-offset] = data[j]
                indices[j-offset] = indices[j]
        else:
            offset += end - start
            offset_rows += 1
    arr.data = data[:nnz]
    arr.indices = indices[:nnz]    
    arr.indptr = indptr[:m-offset_rows]
    arr._shape = (m - offset_rows, n) if m-offset_rows else (0, 0)
    

def filter_sparse_array(arr, ids, metadata, function, axis, invert, remove=True):
    fmt = arr.getformat()
    if fmt not in {'csc', 'csr'}:
        raise TypeError("Format not supported (use CSC/CSR)")

    cdef cnp.ndarray[cnp.int8_t, ndim=1] bools = _make_filter_array(ids, metadata, function, invert)
    
    if axis == 0:
        if fmt == 'csr':
            if remove:
                _remove_rows_csr(arr, bools)
            else:
                _zero_rows_csr(arr, bools)
        elif fmt == 'csc':
            if remove:
                raise TypeError("Compacting array not supported (CSC, axis 0)")
            else:
                _zero_rows_csc(arr, bools)
    elif axis == 1:
        if fmt == 'csr':
            if remove:
                raise TypeError("Compacting array not supported (CSR, axis 1)")
            else:
                _zero_columns_csr(arr, bools)
        elif fmt == 'csc':
            if remove:
                _remove_columns_csc(arr, bools)
            else:
                _zero_columns_csc(arr, bools)
    else:
        raise ValueError("Unsupported axis")

    if remove:
        ids = tuple(compress(ids, bools))
        metadata = tuple(compress(metadata, bools))

    return ids, metadata
