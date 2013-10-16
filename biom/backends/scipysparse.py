#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, BIOM-Format Project"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from itertools import izip
from operator import itemgetter

from numpy import ndarray
from scipy.sparse import coo_matrix

from biom.util import flatten

class ScipySparseMat(object):
    """Based on CSMat implementation by Daniel McDonald."""

    def __init__(self, num_rows, num_cols, dtype=float, data=None):
        # I hate myself for having the empty special case throughout the
        # code... makes it much less elegant. However, scipy doesn't support
        # empty sparse matrices, but we do...
        if num_rows == 0 and num_cols == 0:
            self._matrix = None
        else:
            # TODO: possible optimization is to allow data to be a preexisting
            # scipy.sparse matrix.
            if data is None:
                self._matrix = coo_matrix((num_rows, num_cols), dtype=dtype)
            else:
                # coo_matrix allows zeros to be added as data, and this affects
                # nnz, items, and iteritems. Clean them out here, as this is
                # the only time these zeros can creep in.
                # TODO: do we also want to handle duplicate entries? coo_matrix
                # allows for this, and the entries will be summed when
                # converted, which could be misleading/wrong...
                self._matrix = coo_matrix(data, shape=(num_rows, num_cols),
                                          dtype=dtype)
                self.convert('csr')
                self._matrix.eliminate_zeros()

    def _is_empty(self):
        return self._matrix is None
    is_empty = property(_is_empty)

    def _get_shape(self):
        if self.is_empty:
            return 0,0
        else:
            return self._matrix.shape
    shape = property(_get_shape)

    def _get_dtype(self):
        if self.is_empty:
            return None
        else:
            return self._matrix.dtype
    dtype = property(_get_dtype)

    def _get_format(self):
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
        """

        If ``fmt`` is ``None`` or we're already in the specified format, do
        nothing.
        """
        if not self.is_empty:
            self._matrix = self._matrix.asformat(fmt)

    def transpose(self):
        """Return a transposed copy of self."""
        transposed = self.__class__(self.shape[1], self.shape[0],
                                    dtype=self.dtype)

        if not self.is_empty:
            transposed._matrix = self._matrix.transpose(copy=True)

        return transposed
    T = property(transpose)

    def sum(self):
        return self._matrix.sum()

    def getRow(self, row_idx):
        num_rows, num_cols = self.shape

        if row_idx >= num_rows or row_idx < 0:
            raise IndexError("Row index %d is out of bounds." % row_idx)

        self.convert('csr')

        row_vector = self.__class__(1, num_cols, dtype=self.dtype)
        row_vector._matrix = self._matrix.getrow(row_idx)

        return row_vector

    def getCol(self, col_idx):
        num_rows, num_cols = self.shape

        if col_idx >= num_cols or col_idx < 0:
            raise IndexError("Column index %d is out of bounds." % col_idx)

        self.convert('csc')

        col_vector = self.__class__(num_rows, 1, dtype=self.dtype)
        col_vector._matrix = self._matrix.getcol(col_idx)

        return col_vector

    def items(self):
        """Return [((r,c),v)]. No guaranteed ordering!"""
        return list(self.iteritems())

    def iteritems(self):
        """Generator yielding ((r,c),v). No guaranteed ordering!"""
        if not self.is_empty:
            self.convert('coo')

            for r, c, v in izip(self._matrix.row, self._matrix.col,
                                self._matrix.data):
                yield (r, c), v

    def copy(self):
        """Return a deep copy of self."""
        new_self = self.__class__(*self.shape, dtype=self.dtype)

        if not self.is_empty:
            new_self._matrix = self._matrix.copy()

        return new_self

    def __eq__(self, other):
        """Return True if both matrices are equal.
        
        Matrices are equal iff the following items are equal:
        - type
        - shape
        - dtype
        - size (nnz)
        - matrix data (more expensive, so checked last)

        Sparse format does not need to be the same. ``self`` and ``other`` will
        be converted to csr format if necessary before performing the final
        comparison.
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

            # From http://mail.scipy.org/pipermail/scipy-user/2008-April/016276.html
            # TODO: Do we need abs here?
            if (self._matrix - other._matrix).nnz > 0:
                return False

        return True

    def __ne__(self, other):
        """Return True if both matrices are not equal."""
        return not (self == other)

    def __str__(self):
        if self.is_empty:
            return '<0x0 sparse matrix>'
        else:
            return str(self._matrix)

    def __setitem__(self, args, value):
        if self.is_empty:
            raise IndexError("Cannot set an element on a 0x0 matrix.")

        try:
            row, col = args
        except:
            raise IndexError("Must specify the row and column of the element "
                             "to be set.")

        self.convert('lil')
        if value == 0:
            # TODO: we can support this, but need to watch out for efficiency
            # issues and nnz.
            if self._matrix[row,col] != 0:
                raise ValueError("Cannot set an existing non-zero element to "
                                 "zero.")
        else:
            self._matrix[row,col] = value

    def __getitem__(self, args):
        """Handles slices."""
        if self.is_empty:
            raise IndexError("Cannot retrieve an element from a 0x0 matrix.")

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

            return self._matrix[row,col]


def to_scipy(values, transpose=False, dtype=float):
    """Try to return a populated ScipySparseMat object.

    NOTE: assumes the max value observed in row and col defines the size of the
    matrix.
    """
    # if it is a vector
    if isinstance(values, ndarray) and len(values.shape) == 1:
        if transpose:
            mat = nparray_to_scipy(values[:,newaxis], dtype)
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
    """Convert a list of lists into a ScipySparseMat.

    [[row, col, value], ...]
    """
    rows, cols, values = zip(*data)

    if shape is None:
        n_rows = max(rows) + 1
        n_cols = max(cols) + 1
    else:
        n_rows, n_cols = shape

    # TODO: the CSMat code doesn't respect dtype. Should we pass it here?
    return ScipySparseMat(n_rows, n_cols, data=(values,(rows,cols)))

def nparray_to_scipy(data, dtype=float):
    """Convert a numpy array to a ScipySparseMat."""
    if len(data.shape) == 1:
        shape = (1, data.shape[0])
    else:
        shape = data.shape

    return ScipySparseMat(*shape, dtype=dtype, data=data)

def list_nparray_to_scipy(data, dtype=float):
    """Takes a list of numpy arrays and creates a ScipySparseMat."""
    return ScipySparseMat(len(data), len(data[0]), dtype=dtype, data=data)

def list_scipy_to_scipy(data, dtype=float):
    """Takes a list of ScipySparseMats and creates a ScipySparseMat."""
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
        # TODO: Is this needed here? I thought this was supposed to take a list
        # of sparse matrices, not a list of dicts.
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
    for row_idx,row in enumerate(data):
        # TODO: items() is inefficient. Should use iteritems().
        for (foo,col_idx),val in row.items():
            if is_col:
                # transpose
                rows.append(foo)
                cols.append(row_idx)
                vals.append(val)
            else:
                rows.append(row_idx)
                cols.append(col_idx)
                vals.append(val)

    return ScipySparseMat(n_rows, n_cols, dtype=dtype, data=(vals,(rows,cols)))

# TODO: this seems like the same function as above
def list_dict_to_scipy(data, dtype=float):
    """Takes a list of dict {(row,col):val} and creates a ScipySparseMat."""
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
    for row_idx,row in enumerate(data):
        for (foo,col_idx),val in row.items():
            if is_col:
                # transpose
                rows.append(foo)
                cols.append(row_idx)
                vals.append(val)
            else:
                rows.append(row_idx)
                cols.append(col_idx)
                vals.append(val)

    return ScipySparseMat(n_rows, n_cols, dtype=dtype, data=(vals,(rows,cols)))

def dict_to_scipy(data, dtype=float):
    """Takes a dict {(row,col):val} and creates a ScipySparseMat."""
    # TODO: this shape inferring code should be centralized somewhere. It is
    # heavily duplicated.
    n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
    n_cols = max(data.keys(), key=itemgetter(1))[1] + 1

    rows = []
    cols = []
    vals = []
    for (r,c),v in data.items():
        rows.append(r)
        cols.append(c)
        vals.append(v)

    return ScipySparseMat(n_rows, n_cols, dtype=dtype, data=(vals,(rows,cols)))
