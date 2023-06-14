# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np
cimport numpy as cnp

cdef _subsample_with_replacement(cnp.ndarray[cnp.float64_t, ndim=1] data,
                                 cnp.ndarray[cnp.int32_t, ndim=1] indptr,
                                 cnp.int64_t n,
                                 object rng):
    """Subsample non-zero values of a sparse array with replacement

    Parameters
    ----------
    data : {csr_matrix, csc_matrix}.data
        A 1xM sparse vector data
    indptr : {csr_matrix, csc_matrix}.indptr
        A 1xM sparse vector indptr
    n : int
        Number of items to subsample from `arr`
    rng : Generator instance
        A random generator. This will likely be an instance returned 
        by np.random.default_rng

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
        cnp.int32_t start,end,length
        Py_ssize_t i
        cnp.ndarray[cnp.float64_t, ndim=1] pvals

    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        length = end - start
        counts_sum = data[start:end].sum()
        
        pvals = data[start:end] / counts_sum
        data[start:end] = rng.multinomial(n, pvals)


cdef _subsample_without_replacement(cnp.ndarray[cnp.float64_t, ndim=1] data,
                                    cnp.ndarray[cnp.int32_t, ndim=1] indptr,
                                    cnp.int64_t n,
                                    object rng):
    """Subsample non-zero values of a sparse array w/out replacement

    Parameters
    ----------
    data : {csr_matrix, csc_matrix}.data
        A 1xM sparse vector data
    indptr : {csr_matrix, csc_matrix}.indptr
        A 1xM sparse vector indptr
    n : int
        Number of items to subsample from `arr`
    rng : Generator instance
        A random generator. This will likely be an instance returned 
        by np.random.default_rng

    Returns
    -------
    ndarray
        Subsampled data

    Notes
    -----
    This code was adapted from scikit-bio (`skbio.math._subsample`)

    """
    cdef:
        cnp.int64_t counts_sum, count_el, perm_count_el
        cnp.int64_t count_rem
        cnp.ndarray[cnp.int64_t, ndim=1] permuted
        Py_ssize_t i, idx
        cnp.int32_t length,el,start,end

    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        length = end - start
        counts_sum = data[start:end].sum()
        
        if counts_sum < n:
            data[start:end] = 0
            continue

        permuted = rng.choice(counts_sum, n, replace=False, shuffle=False)
        permuted.sort()

        # now need to do reverse mapping
        # since I am not using np.repeat anymore
        # reminder, old logic was
        #   r = np.arange(length)
        #   unpacked = np.repeat(r, data_i[start:end])
        #   permuted_unpacked = rng.choice(unpacked, n, replace=False, shuffle=False)

        el = 0         # index in result/data
        count_el = 0  # index in permutted
        count_rem = long(data[start])  # since each data has multiple els, sub count there
        data[start] = 0.0
        for idx in range(n):
            perm_count_el = permuted[idx]
            # the array is sorted, so just jump ahead
            while (perm_count_el - count_el) >= count_rem:
               count_el += count_rem
               el += 1
               count_rem = long(data[start+el])
               data[start+el] = 0.0
            count_rem -= (perm_count_el - count_el)
            count_el = perm_count_el

            data[start+el] += 1
        # clean up tail elements
        data[start+el+1:end] = 0.0


def _subsample(arr, n, with_replacement, rng):
    """Subsample non-zero values of a sparse array

    Parameters
    ----------
    arr : {csr_matrix, csc_matrix}
        A 1xM sparse vector
    n : int
        Number of items to subsample from `arr`
    with_replacement : bool
        Whether to permute or use multinomial sampling
    rng : Generator instance
        A random generator. This will likely be an instance returned 
        by np.random.default_rng

    Returns
    -------
    ndarray
        Subsampled data

    Notes
    -----
    This code was adapted from scikit-bio (`skbio.math._subsample`)

    """
    if (with_replacement):
       return _subsample_with_replacement(arr.data, arr.indptr, n, rng)
    else:
       return _subsample_without_replacement(arr.data, arr.indptr, n, rng)
