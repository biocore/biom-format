# __author__ = "Daniel McDonald"
# __copyright__ = "Copyright 2012, BIOM-Format Project"
# __credits__ = ["Daniel McDonald", "Jai Rideout", "Greg Caporaso", 
#                        "Jose Clemente", "Justin Kuczynski"]
# __license__ = "GPL"
# __url__ = "http://biom-format.org"
# __version__ = "1.1.0-dev"
# __maintainer__ = "Daniel McDonald"
# __email__ = "daniel.mcdonald@colorado.edu"

cdef extern from "sparsemat_lib.h" namespace "sparsemat":
    cdef struct items_int:
        # support struct for returning all rows/cols/values from a integer mat
        int *rows
        int *cols
        int *values

    cdef struct items_float:
        # support struct for returning all rows/cols/values from a double 
        # precision floating point mat    
        int *rows
        int *cols
        double *values
    
    cdef cppclass SparseMatFloat:
        # prototype for SparseMatFloat c++ object
        SparseMatFloat()
        void insert(int, int, double)
        double get(int, int)
        void erase(int, int)
        int contains(int, int)
        int length()
        items_float keys()
        items_float items()
        void cleanItems()
        
    cdef cppclass SparseMatInt:
        # prototype for SparseMatInt c++ object
        SparseMatInt()
        void insert(int, int, int)
        int get(int, int)
        void erase(int, int)
        int contains(int, int)
        int length()
        items_int keys()
        items_int items()
        void cleanItems()

# inheritence was not working as expected, so it was scrapped.
cdef class PySparseMatFloat:
    # provide a python interface into the c++ object
    cdef SparseMatFloat *thisptr # so we can refer to self
    cdef int rows
    cdef int cols
    
    def __cinit__(self, int rows, int cols):
        """Using default constructor"""
        if rows < 0:
            raise KeyError("The number of rows must be > 0!")
        if cols < 0:
            raise KeyError("The number of cols must be > 0!")
            
        self.rows = rows
        self.cols = cols
        
        self.thisptr = new SparseMatFloat()

    def __dealloc__(self):
        """Using default destructor"""
        del self.thisptr

    def __str__(self):
        """Output like a python dict"""
        cdef int r
        cdef int c
        cdef double v
        
        output = []
        for (r,c),v in self.items():
            output.append('(%d, %d):%f' % (r,c,v))
        output = ', '.join(output)
        return '{%s}' % output
 
    def _boundcheck(self, int row, int col):
        """Check row/col are sane"""
        if row < 0 or row >= self.rows:
            raise KeyError("Row %d is out of bounds!" % row)
        if col < 0 or col >= self.cols:
            raise KeyError("Col %d is out of bounds!" % col)
             
    def insert(self, int row, int col, double value):
        """Insert a value into the matrix"""
        self._boundcheck(row, col)    
        self.thisptr.insert(<unsigned int>row, <unsigned int>col, value)
   
    def get(self, int row, int col):
        """Get a value from the matrix, return 0.0 if doesn't exist"""
        self._boundcheck(row, col)
        return self.thisptr.get(<unsigned int>row, <unsigned int>col)

    def getRow(self, int row, set row_keys):
        """Bulk get a row, return a new self type"""
        cdef int r
        cdef int c
        cdef double v
        cdef int checked_row
         
        self._boundcheck(row, 0) # assume at least a single column...
        checked_row = <unsigned int> row
        new_row = self.__class__(1, self.cols)
        
        for (r,c) in row_keys:
            new_row.insert(0, c, self.thisptr.get(checked_row, <unsigned int> c))
        
        return new_row
        
    def getCol(self, int col, set col_keys):
        """Bulk get a col, return a new self type"""
        cdef int r
        cdef int c
        cdef double v
        cdef int checked_col
         
        self._boundcheck(0, col) # assume at least a single row...
        checked_col = <unsigned int> col
        new_col = self.__class__(self.rows, 1)
        
        for (r,c) in col_keys:
            new_col.insert(r, 0, self.thisptr.get(<unsigned int> r, checked_col))
        
        return new_col
        
    def erase(self, int row, int col):
        """Delete a value from the matrix"""
        self._boundcheck(row, col)
        self.thisptr.erase(<unsigned int>row, <unsigned int>col)
        
    def contains(self, int row, int col):
        """Return True if the row/col exists, False otherwise"""
        self._boundcheck(row, col)
        if self.thisptr.contains(<unsigned int>row, <unsigned int>col) == 1:
            return True
        else:
            return False
        
    def length(self):
        """Return the number of elements within the matrix"""
        return self.thisptr.length()
        
    def keys(self):
        """Return keys [(row,col)] similar to Python dict.keys()"""
        cdef items_float foo
        cdef Py_ssize_t i
        foo = self.thisptr.keys()
        
        results = []
        for i in range(self.length()):
            results.append((foo.rows[i], foo.cols[i]))
                
        return results
    
    def items(self):
        """Return items [((row,col),value)] similar to Python dict.items()"""
        cdef items_float foo
        cdef Py_ssize_t i
        foo = self.thisptr.items()
        
        results = []
        for i in range(self.length()):
            results.append(((foo.rows[i], foo.cols[i]), foo.values[i]))
                
        return results
    
cdef class PySparseMatInt:
    # provide a python interface into the c++ object
    cdef SparseMatInt *thisptr # so we can refer to self
    cdef int rows
    cdef int cols
    
    def __cinit__(self, int rows, int cols):
        """Using default constructor"""
        if rows < 0:
            raise KeyError("The number of rows must be > 0!")
        if cols < 0:
            raise KeyError("The number of cols must be > 0!")
            
        self.rows = rows
        self.cols = cols
        
        self.thisptr = new SparseMatInt()

    def __dealloc__(self):
        """Using default destructor"""
        del self.thisptr

    def __str__(self):
        """Output like a python dict"""
        cdef int r
        cdef int c
        cdef double v
        output = []
        
        for (r,c),v in self.items():
            output.append('(%d, %d):%d' % (r,c,v))
        output = ', '.join(output)
        return '{%s}' % output
   
    def _boundcheck(self, int row, int col):
        """Check row/col are sane"""
        if row < 0 or row >= self.rows:
            raise KeyError("Row %d is out of bounds!" % row)
        if col < 0 or col >= self.cols:
            raise KeyError("Col %d is out of bounds!" % col)
               
    def insert(self, int row, int col, int value):
        """Insert a value into the matrix"""
        self._boundcheck(row, col)
        self.thisptr.insert(<unsigned int>row, <unsigned int>col, value)
   
    def get(self, int row, int col):
        """Get a value from the matrix, return 0.0 if doesn't exist"""
        self._boundcheck(row, col)
        return self.thisptr.get(<unsigned int>row, <unsigned int>col)

    def getRow(self, int row, set row_keys):
        """Bulk get a row, return a new self type"""
        cdef int r
        cdef int c
        cdef int v
        cdef int checked_row
         
        self._boundcheck(row, 0) # assume at least a single column...
        checked_row = <unsigned int> row
        new_row = self.__class__(1, self.cols)
        
        for (r,c) in row_keys:
            new_row.insert(0, c, self.thisptr.get(checked_row, <unsigned int> c))
        
        return new_row
        
    def getCol(self, int col, set col_keys):
        """Bulk get a col, return a new self type"""
        cdef int r
        cdef int c
        cdef int v
        cdef int checked_col
         
        self._boundcheck(0, col) # assume at least a single row...
        checked_col = <unsigned int> col
        new_col = self.__class__(self.rows, 1)
        
        for (r,c) in col_keys:
            new_col.insert(r, 0, self.thisptr.get(<unsigned int> r, checked_col))
        
        return new_col

    def erase(self, int row, int col):
        """Delete a value from the matrix"""
        self._boundcheck(row, col)
        self.thisptr.erase(<unsigned int>row, <unsigned int>col)

    def contains(self, int row, int col):
        """Return True if the row/col exists, False otherwise"""
        self._boundcheck(row, col)
        if self.thisptr.contains(<unsigned int>row, <unsigned int>col) == 1:
            return True
        else:
            return False
            
    def length(self):
        """Return the number of elements within the matrix"""
        return self.thisptr.length()
        
    def keys(self):
        """Return keys [(row,col)] similar to Python dict.keys()"""
        cdef items_int foo
        cdef Py_ssize_t i
        
        foo = self.thisptr.keys()
        
        results = []
        for i in range(self.length()):
            results.append((foo.rows[i], foo.cols[i]))
                
        return results

    def items(self):
        """Return items [((row,col),value)] similar to Python dict.items()"""
        cdef items_int foo
        cdef Py_ssize_t i
        
        foo = self.thisptr.items()
        
        results = []
        for i in range(self.length()):
            results.append(((foo.rows[i], foo.cols[i]), foo.values[i]))
                
        return results
