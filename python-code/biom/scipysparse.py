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

from scipy.sparse import coo_matrix

class ScipySparseMat(object):
    """Based on CSMat implementation by Daniel McDonald."""

    def __init__(self, num_rows, num_cols, dtype=float, data=None):
        # TODO: possible optimization is to allow data to be a preexisting
        # scipy.sparse matrix.
        if data is None:
            self._matrix = coo_matrix((num_rows, num_cols), dtype=dtype)
        else:
            # TODO: coo_matrix allows zeros to be added as data, and this
            # affects nnz! May want some sanity checks, or make our nnz smarter
            # (e.g., use nonzero() instead, which does seem to work correctly.
            # Or can possibly use eliminate_zeros() or check_format()?
            self._matrix = coo_matrix(data, shape=(num_rows, num_cols),
                                      dtype=dtype)

    def _get_shape(self):
        return self._matrix.shape
    shape = property(_get_shape)

    def _get_dtype(self):
        return self._matrix.dtype
    dtype = property(_get_dtype)

    def _get_format(self):
        return self._matrix.getformat()
    fmt = property(_get_format)

    def _get_size(self):
        """Return the number of non-zero elements (NNZ)."""
        return self._matrix.nnz
    size = property(_get_size)

    def convert(self, fmt=None):
        """If ``fmt`` is ``None`` or we're already in the specified format, do nothing."""
        self._matrix = self._matrix.asformat(fmt)

    def transpose(self):
        """Return a transposed copy of self."""
        transposed = self.__class__(self.shape[1], self.shape[0],
                                    dtype=self.dtype)
        transposed._matrix = self._matrix.transpose(copy=True)
        return transposed
    T = property(transpose)

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
        self.convert('coo')

        for r, c, v in izip(self._matrix.row, self._matrix.col,
                            self._matrix.data):
            yield (r, c), v

    def copy(self):
        """Return a deep copy of self."""
        new_self = self.__class__(*self.shape, dtype=self.dtype)
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
        return str(self._matrix)

    def __setitem__(self, args, value):
        try:
            row, col = args
        except:
            raise IndexError("Must specify the row and column of the element "
                             "to be set.")

        self.convert('csr')
        if value == 0:
            # TODO: we can support this, but need to watch out for efficiency
            # issues and nnz.
            if self._matrix[row,col] != 0:
                raise ValueError("Cannot set an existing non-zero element to "
                                 "zero.")
        else:
            # TODO: may be inefficient for csr/csc (use lil or dok instead?).
            self._matrix[row,col] = value

    def __getitem__(self, args):
        """Handles slices."""
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
            self.convert('csr')
            return self._matrix[row,col]


def to_csmat(values, transpose=False, dtype=float):
    """Tries to returns a populated CSMat object

    NOTE: assumes the max value observed in row and col defines the size of the
    matrix
    """
    # if it is a vector
    if isinstance(values, ndarray) and len(values.shape) == 1:
        if transpose:
            mat = nparray_to_csmat(values[:,newaxis], dtype)
        else:
            mat = nparray_to_csmat(values, dtype)
        return mat
    if isinstance(values, ndarray):
        if transpose:
            mat = nparray_to_csmat(values.T, dtype)
        else:
            mat = nparray_to_csmat(values, dtype)
        return mat
    # the empty list
    elif isinstance(values, list) and len(values) == 0:
        mat = CSMat(0,0)
        return mat
    # list of np vectors
    elif isinstance(values, list) and isinstance(values[0], ndarray):
        mat = list_nparray_to_csmat(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    # list of dicts, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], dict):
        mat = list_dict_to_csmat(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    # list of csmat, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], CSMat):
        mat = list_csmat_to_csmat(values,dtype)
        if transpose:
            mat = mat.T
        return mat
    elif isinstance(values, dict):
        mat = dict_to_csmat(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    elif isinstance(values, CSMat):
        mat = values
        if transpose:
            mat = mat.T
        return mat
    else:
        raise TableException, "Unknown input type"
        
def list_list_to_csmat(data, dtype=float, shape=None):
    """Convert a list of lists into a CSMat

    [[row, col, value], ...]
    """
    rows, cols, values = zip(*data)

    if shape is None:
        n_rows = max(rows) + 1
        n_cols = max(cols) + 1
    else:
        n_rows, n_cols = shape

    mat = CSMat(n_rows, n_cols)
    mat.bulkCOOUpdate(rows, cols, values)
    return mat

def nparray_to_csmat(data, dtype=float):
    """Convert a numpy array to a CSMat"""
    rows = []
    cols = []
    vals = []

    if len(data.shape) == 1:
        mat = CSMat(1, data.shape[0], dtype=dtype)
        for col_idx, val in enumerate(data):
            if val != 0:
                rows.append(0)
                cols.append(col_idx)
                vals.append(val)
    else:
        mat = CSMat(*data.shape, dtype=dtype)
        for row_idx, row in enumerate(data):
            for col_idx, value in enumerate(row):
                if value != 0:
                    rows.append(row_idx)
                    cols.append(col_idx)
                    vals.append(value)
    mat.bulkCOOUpdate(rows, cols, vals)
    return mat

def list_nparray_to_csmat(data, dtype=float):
    """Takes a list of numpy arrays and creates a csmat"""
    mat = CSMat(len(data), len(data[0]),dtype=dtype)
    rows = []
    cols = []
    values = []
    for row_idx, row in enumerate(data):
        if len(row.shape) != 1:
            raise TableException, "Cannot convert non-1d vectors!"
        if len(row) != mat.shape[1]:
            raise TableException, "Row vector isn't the correct length!"

        for col_idx, val in enumerate(row):
            rows.append(row_idx)
            cols.append(col_idx)
            values.append(val)
    mat.bulkCOOUpdate(rows, cols, values)
    return mat

def list_csmat_to_csmat(data, dtype=float):
    """Takes a list of CSMats and creates a CSMat"""
    if isinstance(data[0], CSMat):
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

    mat = CSMat(n_rows, n_cols, dtype=dtype)
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
    mat.bulkCOOUpdate(rows, cols, vals) 
    return mat
    
def list_dict_to_csmat(data, dtype=float):
    """Takes a list of dict {(0,col):val} and creates a CSMat"""
    if isinstance(data[0], CSMat):
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

    mat = CSMat(n_rows, n_cols, dtype=dtype)
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
    
    mat.bulkCOOUpdate(rows, cols, vals)
    return mat

def dict_to_csmat(data, dtype=float):
    """takes a dict {(row,col):val} and creates a CSMat"""
    n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
    n_cols = max(data.keys(), key=itemgetter(1))[1] + 1
    mat = CSMat(n_rows, n_cols,dtype=dtype)
    rows = []
    cols = []
    vals = []

    for (r,c),v in data.items():
        rows.append(r)
        cols.append(c)
        vals.append(v)

    mat.bulkCOOUpdate(rows, cols, vals)
    return mat
