#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, BIOM-Format Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from scipy.sparse import coo_matrix

class ScipySparseMat(object):
    def __init__(self, num_rows, num_cols, dtype=float, data=None):
        # TODO: possible optimization is to allow data to be a preexisting
        # scipy.sparse matrix.
        self.shape = (num_rows, num_cols)
        self.dtype = dtype

        if data is None:
            self._matrix = coo_matrix(self.shape, dtype=self.dtype)
        else:
            # TODO: coo_matrix allows zeros to be added as data, and this
            # affects nnz! May want some sanity checks, or make our nnz smarter
            # (e.g., use nonzero() instead, which does seem to work correctly.
            # Or can possibly use eliminate_zeros() or check_format()?
            self._matrix = coo_matrix(data, shape=self.shape, dtype=self.dtype)

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
        - matrix data (more expensive, so performed last)

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
        """dump priv data"""
        l = []
        l.append(self._order)
        l.append("_coo_values\t" + '\t'.join(map(str, self._coo_values)))
        l.append("_coo_rows\t" + '\t'.join(map(str, self._coo_rows)))
        l.append("_coo_cols\t" + '\t'.join(map(str, self._coo_cols)))
        l.append("_values\t" + '\t'.join(map(str, self._values)))
        l.append("_pkd_ax\t" + '\t'.join(map(str, self._pkd_ax)))
        l.append("_unpkd_ax\t" + '\t'.join(map(str, self._unpkd_ax)))
        return '\n'.join(l)

    def __setitem__(self,args,value):
        """Wrap setitem, complain if out of bounds"""
        try:
            row,col = args
        except:
            # fast support foo[5] = 10, like numpy 1d vectors
            col = args
            row = 0
            args = (row,col)

        if row >= self.shape[0]:
            raise IndexError, "Row %d is out of bounds!" % row
        if col >= self.shape[1]:
            raise IndexError, "Col %d is out of bounds!" % col

        if value == 0:
            if args in self:
                raise ValueError("Cannot set an existing non-zero element to "
                                 "zero.")
        else:
            res = self._getitem(args)
            if res == (None, None, None):
                self._coo_rows.append(row)
                self._coo_cols.append(col)
                self._coo_values.append(value)
            else:
                if self._order == "coo":
                    self._coo_values[res[0]] = value
                else:
                    self._values[res[-1]] = value

    def __getitem__(self,args):
        """Wrap getitem to handle slices"""
        try:
            row,col = args
        except TypeError:
            raise IndexError, "Must specify (row, col)"

        if isinstance(row, slice): 
            if row.start is None and row.stop is None:
                return self.getCol(col)
            else:
                raise AttributeError, "Can only handle full : slices per axis"
        elif isinstance(col, slice):
            if col.start is None and col.stop is None:
                return self.getRow(row)
            else:
                raise AttributeError, "Can only handle full : slices per axis"
        else:
            if row >= self.shape[0] or row < 0:
                raise IndexError, "Row out of bounds!"
            if col >= self.shape[1] or col < 0:
                raise IndexError, "Col out of bounds!"

            res = self._getitem(args)
            if res == (None,None,None):
                return self.dtype(0)
            else:
                if self._order == 'coo':
                    return self._coo_values[res[0]]
                else:
                    return self._values[res[-1]]
                
        return self.dtype(0)

    def _getitem(self, args):
        """Mine for an item
        
        if order is csc | csr, returns
        pkd_ax_idx, unpkd_ax_idx, values_idx 

        if order is coo, returns
        rows_idx, cols_idx, values_idx (all the same thing...)
        """
        if self.hasUpdates():
            self.absorbUpdates()

        row,col = args
        if self._order == 'csr':
            start = self._pkd_ax[row]
            stop = self._pkd_ax[row+1]
            for i,c in enumerate(self._unpkd_ax[start:stop]):
                if c == col:
                    return (row, start+i, start+i)

        elif self._order == 'csc':
            start = self._pkd_ax[col]
            stop = self._pkd_ax[col+1]
            for i,r in enumerate(self._unpkd_ax[start:stop]):
                if r == row:
                    return (start+i, col, start+i)

        elif self._order == "coo":
            # O(N) naive... but likely not a major use case
            idx = 0
            for (r,c) in izip(self._coo_rows, self._coo_cols):
                if r == row and c == col:
                    return (idx, idx, idx)
                idx += 1
        else:
            raise ValueError, "Unknown matrix type: %s" % self._order

        return (None, None, None)

    def _buildCSfromCS(self):
        """Convert csc <-> csr"""
        expanded = self._expand_compressed(self._pkd_ax)
        if self._order == "csr":
            csc = self._toCSC(expanded, self._unpkd_ax, self._values)
            self._pkd_ax, self._unpkd_ax, self._values = csc 
            self._order = "csc"

        elif self._order == "csc":
            csr = self._toCSR(self._unpkd_ax, expanded, self._values)
            self._pkd_ax, self._unpkd_ax, self._values = csr
            self._order = "csr"

    def _expand_compressed(self, pkd_ax):
        """Expands packed axis"""
        expanded = zeros(pkd_ax[-1], dtype=uint32)
        last_idx = 0
        pos = uint32(0)
        for idx in pkd_ax[1:]:
            expanded[last_idx:idx] = pos
            pos += 1
            last_idx = idx
        return expanded
            
    def _buildCOOfromCS(self):
        """Constructs a COO representation from CSC or CSR
        
        Invalidates existing CSC or CSR representation
        """
        coo = self._toCOO(self._pkd_ax,self._unpkd_ax,self._values,self._order)
        coo_rows, coo_cols, coo_values = coo
        self._coo_rows.extend(coo_rows)
        self._coo_cols.extend(coo_cols)
        self._coo_values.extend(coo_values)

        self._values = array([], dtype=self.dtype)
        self._pkd_ax = array([], dtype=uint32)
        self._unpkd_ax = array([], dtype=uint32)
        
        self._order = "coo"

    def _toCOO(self, pkd_ax, unpkd_ax, values, current_order):
        """Returns rows, cols, values"""
        coo_values = list(values)
        expanded_ax = list(self._expand_compressed(pkd_ax))
        
        if current_order == 'csr':
            coo_cols = list(unpkd_ax)
            coo_rows = expanded_ax

        elif current_order == 'csc':
            coo_rows = list(unpkd_ax)
            coo_cols = expanded_ax
        else:
            raise ValueError, "Unknown order: %s" % order

        return (coo_rows, coo_cols, coo_values)
        
    def _buildCSfromCOO(self, order):
        """Build a sparse representation

        order is either csc or csr

        Returns instantly if is stable, throws ValueError if the sparse rep
        is already built
        """
        if order == 'csr':
            csr = self._toCSR(self._coo_rows, self._coo_cols, self._coo_values)
            self._pkd_ax, self._unpkd_ax, self._values = csr
        elif order == 'csc':
            csc = self._toCSC(self._coo_rows, self._coo_cols, self._coo_values)
            self._pkd_ax, self._unpkd_ax, self._values = csc 
        else:
            raise ValueError, "Unknown order: %s" % order

        self._coo_rows = []
        self._coo_cols = []
        self._coo_values = []
        self._order = order

    def _toCSR(self, rows, cols, values):
        """Returns packed_axis, unpacked_axis and values"""
        values = array(values, dtype=self.dtype)
        unpkd_ax = array(cols, dtype=uint32)
        tmp_rows = array(rows, dtype=uint32)

        order = argsort(tmp_rows)
        values = values.take(order)
        tmp_rows = tmp_rows.take(order)
        unpkd_ax = unpkd_ax.take(order)

        v_last = -1
        pkd_ax = []
        pos = 0
        p_last = 0
        expected_row_idx = 0

        # determine starting values idx for each row
        # sort values and columns within each row
        for v in tmp_rows:
            if v != v_last:
                pkd_ax.append(pos)

                # Determine if we skipped any rows (i.e. we have empty rows).
                # Copy the last row's index for all empty rows.
                num_empty_rows = v - expected_row_idx
                if num_empty_rows > 0:
                    pkd_ax.extend([pkd_ax[-1]] * num_empty_rows)
                expected_row_idx = v + 1

                col_order = argsort(unpkd_ax[p_last:pos])
                unpkd_ax[p_last:pos] = unpkd_ax[p_last:pos].take(col_order)
                values[p_last:pos] = values[p_last:pos].take(col_order)
                v_last = v
                p_last = pos
            pos += 1
        # catch last column sort    
        col_order = argsort(unpkd_ax[p_last:pos])
        unpkd_ax[p_last:pos] = unpkd_ax[p_last:pos].take(col_order)
        values[p_last:pos] = values[p_last:pos].take(col_order)
        
        pkd_ax.append(pos)

        num_trailing_rows = self.shape[0] - expected_row_idx
        if num_trailing_rows > 0:
            pkd_ax.extend([pkd_ax[-1]] * num_trailing_rows)

        pkd_ax = array(pkd_ax, dtype=uint32)

        return (pkd_ax, unpkd_ax, values)

    def _toCSC(self, rows, cols, values):
        """Returns packed_axis, unpacked_axis, values"""
        values = array(values, dtype=self.dtype)
        unpkd_ax = array(rows, dtype=uint32)
        tmp_cols = array(cols, dtype=uint32)

        order = argsort(tmp_cols)
        values = values.take(order)
        tmp_cols = tmp_cols.take(order)
        unpkd_ax = unpkd_ax.take(order)

        v_last = -1
        pkd_ax = []
        pos = 0
        p_last = 0
        expected_col_idx = 0

        ### gotta be something in numpy that does this...
        for v in tmp_cols:
            if v != v_last:
                pkd_ax.append(pos)

                num_empty_cols = v - expected_col_idx
                if num_empty_cols > 0:
                    pkd_ax.extend([pkd_ax[-1]] * num_empty_cols)
                expected_col_idx = v + 1

                row_order = argsort(unpkd_ax[p_last:pos])
                unpkd_ax[p_last:pos] = unpkd_ax[p_last:pos].take(row_order)
                values[p_last:pos] = values[p_last:pos].take(row_order)
                v_last = v
                p_last = pos
            pos += 1
        
        # catch last row sort
        row_order = argsort(unpkd_ax[p_last:pos])
        unpkd_ax[p_last:pos] = unpkd_ax[p_last:pos].take(row_order)
        values[p_last:pos] = values[p_last:pos].take(row_order)

        pkd_ax.append(pos)

        num_trailing_cols = self.shape[1] - expected_col_idx
        if num_trailing_cols > 0:
            pkd_ax.extend([pkd_ax[-1]] * num_trailing_cols)

        pkd_ax = array(pkd_ax, dtype=uint32)

        return (pkd_ax, unpkd_ax, values)

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
