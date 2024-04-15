# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


from itertools import compress
from collections.abc import Iterable
from types import FunctionType

import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef cnp.ndarray[cnp.uint8_t, ndim=1] \
    _make_filter_array_general(arr,
                               ids,
                               metadata,
                               func,
                               axis,
                               cnp.uint8_t invert):
    """Faster version of
    [func(vals_i, id_i, md_i) ^ invert for
    (vals_i, id_i, md_i) in zip(ids, metadata, rows/cols)]
    """
    cdef:
        Py_ssize_t i, j, n = arr.shape[::-1][axis]
        cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data, \
                                           row_or_col = np.zeros(n)
        cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr, \
                                         indices = arr.indices
        cnp.ndarray[cnp.uint8_t, ndim=1] bools = \
            np.empty(len(ids), dtype=np.uint8)
        cnp.int32_t start, end

    for i in range(len(ids)):
        start, end = indptr[i], indptr[i+1]

        # The following loop should be equivalent to
        # row_or_col = np.zeros(n)
        # row_or_col.put(indices[start:end], data[start:end])
        for j in range(n):
            if start >= end or j < indices[start]:
                row_or_col[j] = 0
            elif j == indices[start]:
                row_or_col[j] = data[start]
                start += 1

        # After converting the output of the filtering function to a
        # bool, we XOR it with invert (if invert is false it doesn't
        # modify the function output, if it's true it inverts it).
        bools[i] = bool(func(row_or_col, ids[i], metadata[i])) ^ invert

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
    arr._shape = (m - offset_rows, n)

def _filter(arr, ids, metadata, index, ids_to_keep, axis, invert):
    """Filter row/columns of a sparse matrix according to the output of a
    boolean function.

    Parameters
    ----------
    arr : sparse matrix
    ids : 1D array_like
    metadata : 1D array_like
    index : dict
        Maps id to index
    ids_to_keep : function or iterable
    axis : int
    invert : bool

    Returns
    -------
    arr : sparse matrix
    ids : 1D ndarray of dtype object
    metadata : tuple
    """
    invert = bool(invert)
    metadata_is_None = metadata is None

    # General version (i.e., filter functions accepts values, ids and
    # metadata) requires CSR for axis 0 and CSC for axis 1.
    if axis == 0:
        arr = arr.tocsr()
    elif axis == 1:
        arr = arr.tocsc()
    fmt = arr.getformat()

    cdef cnp.ndarray[cnp.uint8_t, ndim=1] bools

    if metadata_is_None:
        metadata = (None,) * len(ids)

    if isinstance(ids_to_keep, Iterable):
        idx = [index[id_] for id_ in ids_to_keep]
        ids_to_keep = np.zeros(len(ids), dtype=bool)
        ids_to_keep.put(idx, True)
        bools = np.bitwise_xor(ids_to_keep, invert).view(np.uint8)
    elif isinstance(ids_to_keep, FunctionType):
        bools = _make_filter_array_general(arr, ids, metadata, ids_to_keep,
                                           axis, invert)
    else:
        raise TypeError("ids_to_keep must be an iterable or a function")

    if axis == 0:
        _remove_rows_csr(arr, bools)
    elif axis == 1:
        arr = arr.T  # arr was CSC, CSR after transposing
        _remove_rows_csr(arr, bools)
        arr = arr.T  # Back to CSC

    ids_dtype = ids.dtype
    ids = np.asarray(list(compress(ids, bools)), dtype=ids_dtype)
    metadata = tuple(compress(metadata, bools))

    if metadata_is_None:
        metadata = None
    return arr, ids, metadata
