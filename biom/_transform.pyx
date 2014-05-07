# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division

import numpy as np
cimport numpy as cnp


def _transform(arr, ids, metadata, function, axis):
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
