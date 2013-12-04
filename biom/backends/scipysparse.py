#!/usr/bin/env python
"""scipy sparse matrix backend"""

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from itertools import izip
from operator import itemgetter

from numpy import asarray, ndarray, newaxis, squeeze
from scipy.sparse import coo_matrix

from biom.exception import TableException
from biom.util import flatten

class ScipySparseMat(object):
    """Sparse matrix backend that utilizes scipy.sparse representations.

    Changes between coo, csr, csc, and lil sparse formats as necessary.

    Based on CSMat implementation by Daniel McDonald.
    """

    @staticmethod
    def convertVectorToDense(vec):
        """Converts a ScipySparseMat row/column vector to a dense numpy array.

        The numpy array that is returned will always be a 1-dimensional row
        vector.
        """
        dense_vec = asarray(vec._matrix.todense())

        if vec.shape == (1, 1):
            # Handle the special case where we only have a single element, but
            # we don't want to return a numpy scalar / 0-d array. We still want
            # to return a vector of length 1.
            return dense_vec.reshape(1)
        else:
            return squeeze(dense_vec)

    def __init__(self, num_rows, num_cols, dtype=float, data=None):
        # I hate myself for having the empty special case throughout the
        # code... makes it much less elegant. However, versions of scipy prior
        # to 0.13.0 don't support empty (i.e., 0x0, 0xn, nx0) sparse matrices,
        # but we do. We'll treat these types of matrices as empty/null but
        # still keep track of their shape.
        if num_rows == 0 or num_cols == 0:
            self._matrix = None
            self._empty_shape = (num_rows, num_cols)
        else:
            if data is None:
                self._matrix = coo_matrix((num_rows, num_cols), dtype=dtype)
            else:
                # coo_matrix allows zeros to be added as data, and this affects
                # nnz, items, and iteritems. Clean them out here, as this is
                # the only time these zeros can creep in.
                # Note: coo_matrix allows duplicate entries; the entries will
                # be summed when converted. Not really sure how we want to
                # handle this generally within BIOM- I'm okay with leaving it
                # as undefined behavior for now.
                self._matrix = coo_matrix(data, shape=(num_rows, num_cols),
                                          dtype=dtype)
                self.convert('csr')
                self._matrix.eliminate_zeros()

    def _is_empty(self):
        """Return ``True`` if the matrix is empty/null.

        A matrix is empty/null if the either the number of rows or columns (or
        both) are zero.
        """
        return self._matrix is None
    is_empty = property(_is_empty)

    def _get_shape(self):
        """Return a two-element tuple indicating the shape of the matrix."""
        if self.is_empty:
            return self._empty_shape
        else:
            return self._matrix.shape
    shape = property(_get_shape)

    def _get_dtype(self):
        """Return the type of data being stored in the matrix."""
        if self.is_empty:
            return None
        else:
            return self._matrix.dtype
    dtype = property(_get_dtype)

    def _get_format(self):
        """Return the current sparse format as a string."""
        if self.is_empty:
            return None
        else:
            return self._matrix.getformat()
    fmt = property(_get_format)

    def _get_size(self):
        """Return the number of non-zero elements (NNZ)."""
        if self.is_empty:
            return 0
        else:
            return self._matrix.nnz
    size = property(_get_size)

    def convert(self, fmt=None):
        """Convert the matrix to the specified sparse format.

        If ``fmt`` is ``None`` or we're already in the specified format, do
        nothing.
        """
        if not self.is_empty:
            self._matrix = self._matrix.asformat(fmt)

    def transpose(self):
        """Return a transposed copy of ``self``."""
        transposed = self.__class__(self.shape[1], self.shape[0],
                                    dtype=self.dtype)

        if not self.is_empty:
            # lil's transpose method doesn't have the copy kwarg, but all of
            # the others do.
            if self.fmt == 'lil':
                self.convert('csr')

            transposed._matrix = self._matrix.transpose(copy=True)

        return transposed
    T = property(transpose)

    def sum(self, axis=None):
        """Sum the entire matrix or along rows/columns.

        ``axis`` can be ``None``, 0, or 1.
        """
        if self.is_empty:
            matrix_sum = 0
        else:
            matrix_sum = squeeze(asarray(self._matrix.sum(axis=axis)))

            # We only want to return a scalar if the whole matrix was summed.
            if axis is not None and matrix_sum.shape == ():
                matrix_sum = matrix_sum.reshape(1)

        return matrix_sum

    def getRow(self, row_idx):
        """Return the row at ``row_idx`` as a ``ScipySparseMat``.

        A row vector will be returned in csr format.
        """
        if self.is_empty:
            raise IndexError("Cannot retrieve a row from an "
                             "empty/null matrix.")

        num_rows, num_cols = self.shape

        if row_idx >= num_rows or row_idx < 0:
            raise IndexError("Row index %d is out of bounds." % row_idx)

        self.convert('csr')

        row_vector = self.__class__(1, num_cols, dtype=self.dtype)
        row_vector._matrix = self._matrix.getrow(row_idx)

        return row_vector

    def getCol(self, col_idx):
        """Return the column at ``col_idx`` as a ``ScipySparseMat``.

        A column vector will be returned in csc format.
        """
        if self.is_empty:
            raise IndexError("Cannot retrieve a column from an empty/null "
                             "matrix.")

        num_rows, num_cols = self.shape

        if col_idx >= num_cols or col_idx < 0:
            raise IndexError("Column index %d is out of bounds." % col_idx)

        self.convert('csc')

        col_vector = self.__class__(num_rows, 1, dtype=self.dtype)
        col_vector._matrix = self._matrix.getcol(col_idx)

        return col_vector

    def items(self):
        """Return ``[((r,c),v)]``. No guaranteed ordering!"""
        return list(self.iteritems())

    def iteritems(self):
        """Generator yielding ``((r,c),v)``. No guaranteed ordering!"""
        if not self.is_empty:
            self.convert('coo')

            for r, c, v in izip(self._matrix.row, self._matrix.col,
                                self._matrix.data):
                yield (r, c), v

    def copy(self):
        """Return a deep copy of ``self``."""
        new_self = self.__class__(*self.shape, dtype=self.dtype)

        if not self.is_empty:
            new_self._matrix = self._matrix.copy()

        return new_self

    def __eq__(self, other):
        """Return ``True`` if both matrices are equal.

        Matrices are equal iff the following items are equal:
        - type
        - shape
        - dtype
        - size (nnz)
        - matrix data (more expensive, so checked last)

        The sparse format does not need to be the same between the two
        matrices. ``self`` and ``other`` will be converted to csr format if
        necessary before performing the final comparison.
        """
        if not isinstance(other, self.__class__):
            return False

        if self.shape != other.shape:
            return False

        if self.dtype != other.dtype:
            return False

        if self.size != other.size:
            return False

        if not self.is_empty:
            self.convert('csr')
            other.convert('csr')

            # From:
            # http://mail.scipy.org/pipermail/scipy-user/2008-April/016276.html
            if abs(self._matrix - other._matrix).nnz > 0:
                return False

        return True

    def __ne__(self, other):
        """Return ``True`` if both matrices are not equal."""
        return not (self == other)

    def __str__(self):
        if self.is_empty:
            return '<%dx%d empty/null sparse matrix>' % self.shape
        else:
            return str(self._matrix)

    def __setitem__(self, args, value):
        """Set the element at the specified row and column.

        Currently does not support setting a non-zero element to zero.
        """
        if self.is_empty:
            raise IndexError("Cannot set an element on an empty/null matrix.")

        try:
            row, col = args
        except:
            raise IndexError("Must specify the row and column of the element "
                             "to be set.")

        self.convert('lil')
        if value == 0:
            # We can support this with scipy.sparse, but need to watch out for
            # efficiency issues and nnz. Leaving this unsupported for now to
            # match CSMat.
            if self._matrix[row, col] != 0:
                raise ValueError("Cannot set an existing non-zero element to "
                                 "zero.")
        else:
            self._matrix[row, col] = value

    def __getitem__(self, args):
        """Handles row or column slices."""
        if self.is_empty:
            raise IndexError("Cannot retrieve an element from an empty/null "
                             "matrix.")

        try:
            row, col = args
        except:
            raise IndexError("Must specify (row, col).")

        if isinstance(row, slice) and isinstance(col, slice):
            raise IndexError("Can only slice a single axis.")

        if isinstance(row, slice):
            if row.start is None and row.stop is None:
                return self.getCol(col)
            else:
                raise IndexError("Can only handle full : slices per axis.")
        elif isinstance(col, slice):
            if col.start is None and col.stop is None:
                return self.getRow(row)
            else:
                raise IndexError("Can only handle full : slices per axis.")
        else:
            if self.fmt == 'coo':
                self.convert('csr')

            return self._matrix[row, col]

def to_scipy(values, transpose=False, dtype=float):
    """Try to return a populated ``ScipySparseMat`` object.

    NOTE: assumes the max value observed in row and col defines the size of the
    matrix.
    """
    # if it is a vector
    if isinstance(values, ndarray) and len(values.shape) == 1:
        if transpose:
            mat = nparray_to_scipy(values[:, newaxis], dtype)
        else:
            mat = nparray_to_scipy(values, dtype)
        return mat
    if isinstance(values, ndarray):
        if transpose:
            mat = nparray_to_scipy(values.T, dtype)
        else:
            mat = nparray_to_scipy(values, dtype)
        return mat
    # the empty list
    elif isinstance(values, list) and len(values) == 0:
        return ScipySparseMat(0, 0)
    # list of np vectors
    elif isinstance(values, list) and isinstance(values[0], ndarray):
        mat = list_nparray_to_scipy(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    # list of dicts, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], dict):
        mat = list_dict_to_scipy(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    # list of ScipySparseMats, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], ScipySparseMat):
        mat = list_scipy_to_scipy(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    elif isinstance(values, dict):
        mat = dict_to_scipy(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    elif isinstance(values, ScipySparseMat):
        mat = values
        if transpose:
            mat = mat.T
        return mat
    else:
        raise TableException("Unknown input type")

def list_list_to_scipy(data, dtype=float, shape=None):
    """Convert a list of lists into a ``ScipySparseMat``.

    [[row, col, value], ...]
    """
    rows, cols, values = izip(*data)

    if shape is None:
        n_rows = max(rows) + 1
        n_cols = max(cols) + 1
    else:
        n_rows, n_cols = shape

    return ScipySparseMat(n_rows, n_cols, data=(values, (rows, cols)))

def nparray_to_scipy(data, dtype=float):
    """Convert a numpy array to a ``ScipySparseMat``."""
    if len(data.shape) == 1:
        shape = (1, data.shape[0])
    else:
        shape = data.shape

    return ScipySparseMat(*shape, dtype=dtype, data=data)

def list_nparray_to_scipy(data, dtype=float):
    """Takes a list of numpy arrays and creates a ``ScipySparseMat``."""
    return ScipySparseMat(len(data), len(data[0]), dtype=dtype, data=data)

def list_scipy_to_scipy(data, dtype=float):
    """Takes a list of ``ScipySparseMat``s and creates a ``ScipySparseMat``."""
    if isinstance(data[0], ScipySparseMat):
        if data[0].shape[0] > data[0].shape[1]:
            is_col = True
            n_cols = len(data)
            n_rows = data[0].shape[0]
        else:
            is_col = False
            n_rows = len(data)
            n_cols = data[0].shape[1]
    else:
        all_keys = flatten([d.keys() for d in data])
        n_rows = max(all_keys, key=itemgetter(0))[0] + 1
        n_cols = max(all_keys, key=itemgetter(1))[1] + 1
        if n_rows > n_cols:
            is_col = True
            n_cols = len(data)
        else:
            is_col = False
            n_rows = len(data)

    rows = []
    cols = []
    vals = []
    for row_idx, row in enumerate(data):
        for (foo, col_idx), val in row.items():
            if is_col:
                # transpose
                rows.append(foo)
                cols.append(row_idx)
                vals.append(val)
            else:
                rows.append(row_idx)
                cols.append(col_idx)
                vals.append(val)

    return ScipySparseMat(n_rows, n_cols, dtype=dtype,
                          data=(vals, (rows, cols)))

def list_dict_to_scipy(data, dtype=float):
    """Takes a list of dict {(row,col):val} and creates a ``ScipySparseMat``."""
    if isinstance(data[0], ScipySparseMat):
        if data[0].shape[0] > data[0].shape[1]:
            is_col = True
            n_cols = len(data)
            n_rows = data[0].shape[0]
        else:
            is_col = False
            n_rows = len(data)
            n_cols = data[0].shape[1]
    else:
        all_keys = flatten([d.keys() for d in data])
        n_rows = max(all_keys, key=itemgetter(0))[0] + 1
        n_cols = max(all_keys, key=itemgetter(1))[1] + 1
        if n_rows > n_cols:
            is_col = True
            n_cols = len(data)
        else:
            is_col = False
            n_rows = len(data)

    rows = []
    cols = []
    vals = []
    for row_idx, row in enumerate(data):
        for (foo, col_idx), val in row.items():
            if is_col:
                # transpose
                rows.append(foo)
                cols.append(row_idx)
                vals.append(val)
            else:
                rows.append(row_idx)
                cols.append(col_idx)
                vals.append(val)

    return ScipySparseMat(n_rows, n_cols, dtype=dtype,
                          data=(vals, (rows, cols)))

def dict_to_scipy(data, dtype=float):
    """Takes a dict {(row,col):val} and creates a ``ScipySparseMat``."""
    n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
    n_cols = max(data.keys(), key=itemgetter(1))[1] + 1

    rows = []
    cols = []
    vals = []
    for (r, c), v in data.items():
        rows.append(r)
        cols.append(c)
        vals.append(v)

    return ScipySparseMat(n_rows, n_cols, dtype=dtype,
                          data=(vals, (rows, cols)))
