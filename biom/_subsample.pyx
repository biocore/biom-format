# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# distutils: language = c++
#
import numpy as np
cimport numpy as cnp

cdef extern from "_subsample_cpp.cpp":
    pass

cdef extern from "_subsample_cpp.hpp":
    cdef cppclass WeightedSample:
        WeightedSample(unsigned int _max_count, unsigned int _n, unsigned int random_seed)
        void do_sample(double* data_arr, int start, int end)

def _subsample_with_replacement(arr, n, rng):
    """Subsample non-zero values of a sparse array with replacement

    Parameters
    ----------
    arr : {csr_matrix, csc_matrix}
        A 1xM sparse vector
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
        cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
        cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
        Py_ssize_t i, length

    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        length = end - start
        counts_sum = data[start:end].sum()
        
        pvals = data[start:end] / counts_sum
        data[start:end] = rng.multinomial(n, pvals)


def _subsample_without_replacement(arr, n, rng):
    """Subsample non-zero values of a sparse array w/out replacement

    Parameters
    ----------
    arr : {csr_matrix, csc_matrix}
        A 1xM sparse vector
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
        cnp.ndarray[cnp.float64_t, ndim=1] data = arr.data
        cnp.ndarray[cnp.int32_t, ndim=1] indptr = arr.indptr
        cnp.ndarray[cnp.int32_t, ndim=1] lengths
        Py_ssize_t i
        cnp.uint32_t length,max_len
        cnp.uint32_t cn = n
        WeightedSample *sample_data

    lengths = np.empty(indptr.shape[0] - 1, dtype=np.int32)
    for i in range(indptr.shape[0] - 1):
        start, end = indptr[i], indptr[i+1]
        length = end - start
        lengths[i] = length
        counts_sum = data[start:end].sum()
        if counts_sum < n:
           data[start:end] = 0
           length = 0 # special value to signal to skip
        lengths[i] = length

    max_len = lengths.max()
    
    sample_data = new WeightedSample(max_len, cn,
                                     rng.integers(0,2**32, dtype=np.uint32))
    for i in range(indptr.shape[0] - 1):
        if lengths[i]==0:
            continue
        sample_data.do_sample(&data[0], indptr[i], indptr[i+1])


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
       return _subsample_with_replacement(arr, n, rng)
    else:
       return _subsample_without_replacement(arr, n, rng)
