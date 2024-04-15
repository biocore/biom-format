# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


import numpy as np
cimport numpy as cnp
cnp.import_array()


def _transform(arr, ids, metadata, function, axis):
    """Transform non-zero values of a sparse array, in place.

    Only non null values can be modified: the density of the sparse
    matrix can't increase. However, zeroing values is fine.

    Parameters
    ----------
    arr : csr_matrix or csc_matrix
        Matrix whose rows or columns (respectively) are to be
        transformed.
    ids : 1D array_like
        ids along the given axis.
    metadata : 1D array_like or None
        metadata along the given axis.
    function : function
        A function that takes three values: an array of nonzero values
        for each column or row of `arr`, an id string and a metadata
        dictionary. It must return an array of transformed values.
    axis : int
        Transform rows of `arr` if 0, columns if 1.
    """
    cdef:
        Py_ssize_t n, row_or_col
        cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
        cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
        cdef cnp.int32_t start, end

    if metadata is None:
        metadata = (None,) * len(ids)

    n = arr.shape[axis]
    for row_or_col in range(n):
        start, end = indptr[row_or_col], indptr[row_or_col+1]
        id_ = ids[row_or_col]
        md = metadata[row_or_col]
        data[start:end] = function(data[start:end], id_, md)
