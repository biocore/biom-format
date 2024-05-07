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


cdef _subsample_with_replacement(cnp.ndarray[cnp.float64_t, ndim=1] data,
                                 cnp.ndarray[cnp.int32_t, ndim=1] indptr,
                                 cnp.int64_t n,
                                 object rng):
    """Subsample non-zero values of a sparse array with replacement

    Note: this method operates in place

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

    Notes
    -----
    This code was adapted from scikit-bio (`skbio.math._subsample`)

    """
    cdef:
        cnp.float64_t counts_sum
        cnp.int32_t start,end,length
        Py_ssize_t i
        cnp.ndarray[cnp.float64_t, ndim=1] pvals
        cnp.ndarray[cnp.float64_t, ndim=1] data_ceil 
        
    data_ceil = np.ceil(data)
    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        length = end - start

        # base p-values on integer data to avoid small numerical issues with 
        # float on sum
        counts_sum = data_ceil[start:end].sum()
        pvals = data_ceil[start:end] / counts_sum

        data[start:end] = rng.multinomial(n, pvals)


cdef _subsample_without_replacement(cnp.ndarray[cnp.float64_t, ndim=1] data,
                                    cnp.ndarray[cnp.int32_t, ndim=1] indptr,
                                    cnp.int64_t n,
                                    object rng):
    """Subsample non-zero values of a sparse array w/out replacement

    Note: this method operates in place

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
    """
    cdef:
        cnp.int64_t counts_sum, count_el, perm_count_el
        cnp.int64_t count_rem
        cnp.ndarray[cnp.int64_t, ndim=1] permuted, intdata
        Py_ssize_t i, idx
        cnp.int32_t length,el,start,end
        cnp.int64_t el_cnt

    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        length = end - start
        # We are relying on data being integers
        # If there are rounding erros, fp64 sums can lead to
        # big errors in sum, so convert to int64, first
        intdata = data[start:end].astype(np.int64)
        counts_sum = intdata.sum()
        
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
        # 
        # specifically, what we're going to do here is randomly pick what elements within
        # each sample to keep. this is analogous issuing the prior np.repeat call, and obtaining
        # a random set of index positions for that resulting array. however, we do not need to 
        # perform the np.repeat call as we know the length of that resulting vector already,
        # and additionally, we can compute the sample associated with an index in that array
        # without constructing it.

        el = 0         # index in result/data
        count_el = 0  # index in permutted
        count_rem = intdata[0]  # since each data has multiple els, keep track how many are left
        el_cnt = 0
        for idx in range(n):
            perm_count_el = permuted[idx]
            # The array is sorted, so just jump ahead if needed
            # Move until we get withing the elements range
            while (perm_count_el - count_el) >= count_rem:
               #save the computed value
               data[start+el] = el_cnt
               # move to next element
               el += 1
               # move to the beginning of next element
               count_el += count_rem
               # Load how much we have avaialble
               count_rem = intdata[el]
               #re-start the el counter
               el_cnt = 0
            # increment the el counter
            el_cnt += 1
            # update the counters
            # reduce what is left
            count_rem -= (perm_count_el - count_el)
            #move the pointer to where we stopped
            count_el = perm_count_el
        # save the last value
        data[start+el] = el_cnt
        # clean up tail elements
        data[start+el+1:end] = 0


def subsample(arr, n, with_replacement, rng):
    """Subsample non-zero values of a sparse array

    Note: this method operates in place

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

    Notes
    -----
    This code was adapted from scikit-bio (`skbio.math._subsample`)

    """
    if (with_replacement):
       _subsample_with_replacement(arr.data, arr.indptr, n, rng)
    else:
       _subsample_without_replacement(arr.data, arr.indptr, n, rng)
