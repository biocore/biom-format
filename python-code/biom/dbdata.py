#!/usr/bin/env python

from biomdb import BiomDB
from numpy import zeros
from biom.table import Table
biomdb = BiomDB()

class dbData(object):
    """ """
    _DB_TYPES = {int:'INT',float:'FLOAT'}
    def __init__(self, observations, samples, dtype=float):
        self.Observations = observations
        self.Samples = samples
        self.shape = (len(observations), len(samples))
        self.dtype = self._DB_TYPES[dtype] # casting is minimal, trust the programmer...
        self.Table = biomdb.createTempTable(dtype=self.dtype)
        
        self._obs_r_index = dict([(i,o) for i,o in enumerate(observations)])
        self._obs_index = dict([(o,i) for i,o in enumerate(observations)])
        self._sample_r_index = dict([(i,s) for i,s in enumerate(samples)])
        self._sample_index = dict([(s,i) for i,s in enumerate(samples)])
        
        assert len(self._obs_index) == self.shape[0], "Observation ids might not be unique!"
        assert len(self._sample_index) == self.shape[1], "Sample ids might not be unique!"
    
    def __del__(self):
        """Toss the table associated with this object apon deletion"""
        #biomdb.dropTable(self.Table)
        print "not dropping table %s" % self.Table
        pass
    def __eq__(self, other):
        """checks for eq"""
        return biomdb.tableEquality(self.Table, other.Table)
            
    def __setitem__(self,args,value):
        """Wrap setitem, complain if out of bounds"""
        row,col = args
        in_self_rows, in_self_cols = self.shape

        if row >= in_self_rows or row < 0:
            raise IndexError, "The specified row is out of bounds"
        if col >= in_self_cols or col < 0:
            raise IndexError, "The specified col is out of bounds"

        biomdb.setItem(self.Table, self._obs_r_index[row], 
                                 self._sample_r_index[col], value)

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
            # boundary check
            self_rows, self_cols = self.shape
            if row >= self_rows or row < 0:
                raise IndexError, "Row index out of range"
            if col >= self_cols or col < 0:
                raise IndexError, "Col index out of range"

            return biomdb.getItem(self.Table, self._obs_r_index[row], 
                                         self._sample_r_index[col])
            
    def getRow(self, row):
        """Returns a row: {((row,col):value}"""
        in_self_rows, in_self_cols = self.shape
        if row >= in_self_rows or row < 0:
            raise IndexError, "The specified row is out of bounds"
    
        row_slice = zeros((1,in_self_cols))
        
        for r,c,v in biomdb.getSampleByObs(self.Table, self._obs_r_index[row]):
            row_slice[0, self._sample_index[c]] = v
        
        return row_slice
        
    def getCol(self, col):
        """Return a col: {((row,col):value}"""
        in_self_rows, in_self_cols = self.shape
        if col >= in_self_cols or col < 0:
            raise IndexError, "The specified col is out of bounds"

        col_slice = zeros((1, in_self_rows))
        
        for r,c,v in biomdb.getObsBySample(self.Table, self._sample_r_index[col]):
            col_slice[0, self._obs_index[r]]
        
        return col_slice
        
    def copy(self):
        """Copy self"""
        new_self = self.__class__(self.Observations, self.Samples)
        new_table = biomdb.createTempTable(from_table=self.Table, dtype=self.dtype)
        new_self.Table = new_table
        return new_self

    def transpose(self):
        """Transpose self"""
        print "in transpose"
        print self.getRow(0)
        print self.getRow(1)
        print self.Table
        new_self = self.__class__(self.Samples, self.Observations)
        new_table = biomdb.createTempTable(from_table=self.Table, dtype=self.dtype)
        new_self.Table = new_table
        print new_self.getRow(0)
        print new_self.getRow(1)
        
        return new_self
        
class dbTable(Table):
    _biom_matrix_type = "db"
    
    def __init__(*args, **kwargs):
        super(dbTable, self).__init__(*args,**kwargs)
        
    def _data_equality(self, other):
        raise NotImplementedError
    
    def _conv_to_np(self, v):
        raise NotImplementedError
        
    def conv_to_self_type(self, vals, transpose=False):
        raise NotImplementedError
    
    def __iter__(self):
        return self.iterSamples()
        
    def _iter_samp(self):
        raise NotImplementedError
    
    def _iter_obs(self):
        raise NotImplementedError
