# ----------------------------------------------------------------------------
# Copyright (c) 2013--, biom development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

import numpy as np
cimport numpy as cnp


def _subsample(arr, n):
    """Subsample non-zero values of a sparse array

    Parameters
    ----------
    arr : {csr_matrix, csc_matrix}
        A 1xM sparse vector
    n : int
        Number of items to subsample from `arr`
    
    Returns
    -------
    ndarray
        Subsampled data

    Notes
    -----
    This code was adapted from scikit-bio (`skbio.math._subsample`)

    """
    cdef:
        cnp.int64_t counts_sum
        cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
        cnp.ndarray[cnp.float64_t, ndim=1] result
        cnp.ndarray[cnp.int32_t, ndim=1] indices = arr.indices
        cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
        cnp.ndarray[cnp.int32_t, ndim=1] permuted, unpacked
        cnp.float64_t cnt
        Py_ssize_t unpacked_idx, i, j

    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        counts_sum = data[start:end].sum()
       
        if counts_sum < n:
            data[start:end] = 0
            continue

        unpacked = np.empty(counts_sum, dtype=np.int32)
        unpacked_idx = 0

        for i in range(start, end):
            cnt = data[i]

            for j in range(int(cnt)):
                unpacked[unpacked_idx] = i - start
                unpacked_idx += 1
       
        permuted = np.random.permutation(unpacked)[:n]
        
        result = np.zeros(end - start, dtype=np.float64)
        for idx in range(permuted.shape[0]):
            result[permuted[idx]] += 1

        data[start:end] = result
