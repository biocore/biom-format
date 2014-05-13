# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division

from itertools import compress
from collections import Iterable
from types import FunctionType

import numpy as np
cimport numpy as cnp

from biom.exception import TableException


cdef _zero_rows_CSR_or_columns_CSC(arr,
                                   cnp.ndarray[cnp.uint8_t, ndim=1] booleans,
                                   int axis):
    """Zero out rows or columns for a matrix in CSR or CSC format
    respectively."""
    cdef Py_ssize_t n, row_or_col
    cdef cnp.int32_t start, end
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    n = arr.shape[axis]
    for row_or_col in range(n):
        if booleans[row_or_col]:
            continue
        start, end = indptr[row_or_col], indptr[row_or_col+1]
        data[start:end] = 0
    arr.eliminate_zeros()

cdef _zero_columns_CSR_or_rows_CSC(arr,
                                   cnp.ndarray[cnp.uint8_t, ndim=1] booleans):
    """Zero out rows or columns for a matrix in CSC or CSR format
    respectively."""
    cdef Py_ssize_t i, col_or_row
    cdef cnp.ndarray[cnp.int32_t, ndim=1] indices = arr.indices
    cdef cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
    for i in range(arr.indices.size):
        col_or_row = indices[i]
        if booleans[col_or_row]:
            continue
        data[i] = 0
    arr.eliminate_zeros()

cdef cnp.ndarray[cnp.uint8_t, ndim=1] \
    _make_filter_array_general(arr,
                               ids,
                               metadata,
                               func,
                               cnp.uint8_t invert):
    """Faster version of
    [func(id_i, md_i, vals_i) ^ invert for
    (id_i, md_i, vals_i) in zip(ids, metadata, rows/cols)]
    """
    cdef:
        cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
        cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
        cnp.ndarray[cnp.uint8_t, ndim=1] bools = \
            np.empty(len(ids), dtype=np.uint8)
        cnp.int32_t start, end
        Py_ssize_t i
    for i in range(len(ids)):
        start, end = indptr[i], indptr[i+1]
        bools[i] = bool(func(ids[i], metadata[i], data[start:end])) ^ invert
    return bools

cdef _remove_rows_csr(arr, cnp.ndarray[cnp.uint8_t, ndim=1] booleans):
    """Sparse equivalent of arr[booleans] for a dense array.
    """
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
    arr.indptr = indptr[:m-offset_rows+1]
    arr._shape = (m - offset_rows, n) if m-offset_rows else (0, 0)

def _filter(arr, ids, metadata, ids_to_keep, axis, invert,
            remove=True):
    """Filter row/columns of a sparse matrix according to the output of a
    boolean function.

    Parameters
    ----------
    arr : sparse matrix
    ids : 1D array_like
    metadata : 1D array_like
    ids_to_keep : function or iterable
    axis : int
    invert : bool
    remove : bool
        Whether to "compact" or not the filtered matrix (i.e., keep
        the original size if ``False``, else reduce the shape of the
        returned matrix.

    Returns
    -------
    arr : sparse matrix
    ids : 1D ndarray of dtype object
    metadata : tuple
    """
    invert = bool(invert)
    metadata_is_None = metadata is None

    # General version (i.e., filter functions accepts ids, metadata
    # and values) requires CSR for axis 0 and CSC for axis 1.
    if axis == 0:
        arr = arr.tocsr()
    elif axis == 1:
        arr = arr.tocsc()
    fmt = arr.getformat()

    cdef cnp.ndarray[cnp.uint8_t, ndim=1] bools

    if isinstance(ids_to_keep, Iterable):
        bools = np.bitwise_xor(np.asarray(ids_to_keep, dtype=bool),
                               invert).view(np.uint8)
    elif isinstance(ids_to_keep, FunctionType):
        if metadata_is_None:
            metadata = (None,) * len(ids)

        bools = _make_filter_array_general(arr, ids, metadata, ids_to_keep,
                                           invert)
    else:
        raise TypeError("ids_to_keep must be an iterable or a function")

    if np.all(bools == 0):
        raise TableException("All data was filtered out!")

    if axis == 0:
        if remove:
            _remove_rows_csr(arr, bools)
        else:
            if fmt == 'csr':
                 _zero_rows_CSR_or_columns_CSC(arr, bools, axis)
            elif fmt == 'csc':
                _zero_columns_CSR_or_rows_CSC(arr, bools)
    elif axis == 1:
        if remove:
            arr = arr.T  # arr was CSC, CSR after transposing
            _remove_rows_csr(arr, bools)
            arr = arr.T  # Back to CSC
        else:
            if fmt == 'csr':
                _zero_columns_CSR_or_rows_CSC(arr, bools)
            elif fmt == 'csc':
                 _zero_rows_CSR_or_columns_CSC(arr, bools, axis)

    if remove:
        ids = np.asarray(list(compress(ids, bools)), dtype=object)
        metadata = tuple(compress(metadata, bools))

    if metadata_is_None:
        metadata = None
    return arr, ids, metadata
