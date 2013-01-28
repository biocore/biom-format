#!/usr/bin/env python

from numpy import array, ndarray
from operator import itemgetter
from biom.util import flatten
from biom.exception import TableException
import _sparsemat

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso", 
               "Jose Clemente", "Justin Kuczynski"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.2"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class SparseMat():
    """Support wrapper for the light weight c++ sparse mat

    Must specify rows and columns in advance

    Object cannot "grow" in shape

    There is additional overhead on inserts in order to support rapid lookups
    across rows or columns
    """

    def __init__(self, rows, cols, dtype=float, enable_indices=True):
        self.shape = (rows, cols) 
        self.dtype = dtype # casting is minimal, trust the programmer...

        if dtype == float:
            self._data = _sparsemat.PySparseMatFloat(rows, cols)
        elif dtype == int:
            self._data = _sparsemat.PySparseMatInt(rows, cols)
        else:
            raise ValueError, "Unsupported type %s for _sparsemat" % dtype
            
        if enable_indices:
            self._index_rows = [set() for i in range(rows)]
            self._index_cols = [set() for i in range(cols)]
        else:
            self._index_rows = None
            self._index_cols = None
        self._indices_enabled = enable_indices

    def rebuildIndices(self):
        """(re)Build indices"""
        ir = [set() for i in range(self.shape[0])]
        ic = [set() for i in range(self.shape[1])]

        for ((r,c),v) in self.items():
            ir[r].add((r,c))
            ic[c].add((r,c))

        self._index_rows = ir
        self._index_cols = ic

    def items(self):
        """Generater returning ((r,c),v)"""
        return self._data.items()
    
    def iteritems(self):
        for x in self._data.items():
            yield x

    def __contains__(self, args):
        """Return True if args are in self, false otherwise"""
        row, col = args
        if self._data.contains(row, col):
            return True
        else:
            return False
    
    def erase(self, row, col):
        """Deletes the item at row,col"""
        self._update_internal_indices((row,col), 0)
        self._data.erase(row, col)
             
    def copy(self):
        """Return a copy of self"""
        new_self = self.__class__(self.shape[0], self.shape[1], self.dtype, \
                                  self._indices_enabled)
        for k,v in self._data.items():
            new_self[k] = v
        return new_self

    def __eq__(self, other):
        """Returns true if both SparseMats are the same"""
        if self.shape != other.shape:
            return False
            
        self_keys = set(self._data.keys())
        other_keys = set(other._data.keys())
        
        if len(self_keys) == 0 and len(other_keys) == 0:
            return True
            
        if self_keys != other_keys:
            return False
        
        for k in self_keys:
            if self[k] != other[k]:
                return False
        
        return True
    
    def __ne__(self, other):
        """Return true if both SparseMats are not equal"""
        return not (self == other)
           
    def __setitem__(self,args,value):
        """Wrap setitem, complain if out of bounds"""
        try:
            row,col = args
        except:
            # fast support foo[5] = 10, like numpy 1d vectors
            col = args
            row = 0
            args = (row,col) # passed onto update_internal_indices
            
        if value == 0:
            if args in self:
                self.erase(row, col)
            else:
                return
        else:
            self._update_internal_indices(args, value)
            self._data.insert(row,col,value)
            
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
            return self._data.get(row,col)

    def _update_internal_indices(self, args, value):
        """Update internal row,col indices"""
        if not self._indices_enabled:
            return
        
        row,col = args
        if value == 0:
            if args in self:
                try:
                    self._index_rows[row].remove(args)
                    self._index_cols[col].remove(args)
                except IndexError:
                    raise KeyError, "Row or col is out of bounds!"
            else:
                return # short circuit, no point in setting 0
        else:
            try:
                self._index_rows[row].add(args)
                self._index_cols[col].add(args)
            except IndexError:
                raise KeyError, "Row or col is out of bounds!"

    def getRow(self, row):
        """Returns a row in SparseMat form"""
        in_self_rows, in_self_cols = self.shape
        new_row = SparseMat(1, in_self_cols, enable_indices=False, \
                            dtype=self.dtype)
        if row >= in_self_rows:
            raise IndexError, "Row is out of bounds"

        new_row._data = self._data.getRow(row, self._index_rows[row])
        
        return new_row

    def getCol(self, col):
        """Return a col in SparseMat form"""
        in_self_rows, in_self_cols = self.shape
        new_col = SparseMat(in_self_rows, 1, enable_indices=False, \
                            dtype=self.dtype)
        if col >= in_self_cols:
            raise IndexError, "Col is out of bounds"
        new_col._data = self._data.getCol(col, self._index_cols[col])
        
        return new_col

    def transpose(self):
        """Transpose self"""
        new_self = self.__class__(*self.shape[::-1], \
                enable_indices=self._indices_enabled, dtype=self.dtype)
        for (r,c),v in self._data.items():
            new_self[c,r] = v
        return new_self
    T = property(transpose)

    def _get_size(self):
        """Returns the number of nonzero elements stored"""
        return self._data.length()
    size = property(_get_size)

    def update(self, update_dict):
        """Update self"""
        in_self_rows, in_self_cols = self.shape
    
        for (row,col),value in update_dict.items():
            self[row,col] = value
            
            if self._indices_enabled:
                self._update_internal_indices((row,col), value)

def to_sparsemat(values, transpose=False, dtype=float):
    """Tries to returns a populated SparseMat object

    NOTE: assumes the max value observed in row and col defines the size of the
    matrix
    """
    # if it is a vector
    if isinstance(values, ndarray) and len(values.shape) == 1:
        if transpose:
            mat = nparray_to_sparsemat(values[:,newaxis], dtype)
        else:
            mat = nparray_to_sparsemat(values, dtype)
        return mat
    if isinstance(values, ndarray):
        if transpose:
            mat = nparray_to_sparsemat(values.T, dtype)
        else:
            mat = nparray_to_sparsemat(values, dtype)
        return mat 
    # the empty list
    elif isinstance(values, list) and len(values) == 0:
        mat = SparseMat(0,0)
        return mat
    # list of np vectors
    elif isinstance(values, list) and isinstance(values[0], ndarray):
        mat = list_nparray_to_sparsemat(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    # list of dicts, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], dict):
        mat = list_dict_to_sparsemat(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    # list of sparsemat, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], SparseMat):
        mat = list_sparsemat_to_sparsemat(values,dtype)
        if transpose:
            mat = mat.T
        return mat
    elif isinstance(values, dict):
        mat = dict_to_sparsemat(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    else:
        raise TableException, "Unknown input type"
        
def list_list_to_sparsemat(data, dtype=float, shape=None):
    """Convert a list of lists into a SparseMat

    [[row, col, value], ...]
    """
    d = dict([((r,c),v) for r,c,v in data if v != 0])

    if shape is None:
        n_rows = 0
        n_cols = 0
        for (r,c) in d:
            if r >= n_rows:
                n_rows = r + 1 # deal with 0-based indexes
            if c >= n_cols:
                n_cols = c + 1

        mat = SparseMat(n_rows, n_cols, dtype=dtype)
    else:
        mat = SparseMat(*shape, dtype=dtype)

    mat.update(d)
    return mat

def nparray_to_sparsemat(data, dtype=float):
    """Convert a numpy array to a SparseMat"""
    if len(data.shape) == 1:
        mat = SparseMat(1, data.shape[0], dtype=dtype, enable_indices=False)

        for idx,v in enumerate(data):
            if v != 0:
                mat[(0,idx)] = v
    else:
        mat = SparseMat(*data.shape)
        for row_idx, row in enumerate(data):
            for col_idx, value in enumerate(row):
                if value != 0:
                    mat[(row_idx, col_idx)] = value
    return mat

def list_nparray_to_sparsemat(data, dtype=float):
    """Takes a list of numpy arrays and creates a SparseMat"""
    mat = SparseMat(len(data), len(data[0]),dtype=dtype)
    for row_idx, row in enumerate(data):
        if len(row.shape) != 1:
            raise TableException, "Cannot convert non-1d vectors!"
        if len(row) != mat.shape[1]:
            raise TableException, "Row vector isn't the correct length!"

        for col_idx, val in enumerate(row):
            mat[row_idx, col_idx] = val
    return mat

def list_sparsemat_to_sparsemat(data, dtype=float):
    """Takes a list of SparseMats and creates a SparseMat"""
    if isinstance(data[0], SparseMat):
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

    mat = SparseMat(n_rows, n_cols,dtype=dtype)
    for row_idx,row in enumerate(data):
        for (foo,col_idx),val in row.items():
            if is_col:
                mat[foo,row_idx] = val
            else:
                mat[row_idx,col_idx] = val

    return mat
    
def list_dict_to_sparsemat(data, dtype=float):
    """Takes a list of dict {(0,col):val} and creates a SparseMat"""
    if isinstance(data[0], SparseMat):
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

    mat = SparseMat(n_rows, n_cols, dtype=dtype)
    for row_idx,row in enumerate(data):
        for (foo,col_idx),val in row.items():
            if is_col:
                mat[foo,row_idx] = val
            else:
                mat[row_idx,col_idx] = val

    return mat

def dict_to_sparsemat(data, dtype=float):
    """takes a dict {(row,col):val} and creates a SparseMat"""
    n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
    n_cols = max(data.keys(), key=itemgetter(1))[1] + 1
    mat = SparseMat(n_rows, n_cols,dtype=dtype)
    mat.update(data)
    return mat
