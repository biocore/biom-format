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

# inheritence was not working as expected, so it was scrapped.
cdef class PySparseMatFloat:
    # provide a python interface into the c++ object
    cdef SparseMatFloat *thisptr # so we can refer to self
    
    def __cinit__(self):
        """Using default constructor"""
        self.thisptr = new SparseMatFloat()

    def __dealloc__(self):
        """Using default destructor"""
        del self.thisptr

    def __str__(self):
        """Output like a python dict"""
        output = []
        for (r,c),v in self.items():
            output.append('(%d, %d):%f' % (r,c,v))
        output = ', '.join(output)
        return '{%s}' % output
        
    def insert(self, row, col, value):
        """Insert a value into the matrix"""
        self.thisptr.insert(<unsigned int>row, <unsigned int>col, value)
   
    def get(self, row, col):
        """Get a value from the matrix, return 0.0 if doesn't exist"""
        return self.thisptr.get(<unsigned int>row, <unsigned int>col)

    def erase(self, row, col):
        """Delete a value from the matrix"""
        self.thisptr.erase(<unsigned int>row, <unsigned int>col)
        
    def contains(self, row, col):
        """Return True if the row/col exists, False otherwise"""
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
    
    def __cinit__(self):
        """Using default constructor"""
        self.thisptr = new SparseMatInt()

    def __dealloc__(self):
        """Using default destructor"""
        del self.thisptr

    def __str__(self):
        """Output like a python dict"""
        output = []
        for (r,c),v in self.items():
            output.append('(%d, %d):%d' % (r,c,v))
        output = ', '.join(output)
        return '{%s}' % output
   
    def insert(self, row, col, value):
        """Insert a value into the matrix"""
        self.thisptr.insert(<unsigned int>row, <unsigned int>col, value)
   
    def get(self, row, col):
        """Get a value from the matrix, return 0.0 if doesn't exist"""
        return self.thisptr.get(<unsigned int>row, <unsigned int>col)

    def erase(self, row, col):
        """Delete a value from the matrix"""
        self.thisptr.erase(<unsigned int>row, <unsigned int>col)

    def contains(self, row, col):
        """Return True if the row/col exists, False otherwise"""
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
        foo = self.thisptr.keys()
        
        results = []
        for i in range(self.length()):
            results.append((foo.rows[i], foo.cols[i]))
        
        return results

    def items(self):
        """Return items [((row,col),value)] similar to Python dict.items()"""
        cdef items_int foo
        foo = self.thisptr.items()
        
        results = []
        for i in range(self.length()):
            results.append(((foo.rows[i], foo.cols[i]), foo.values[i]))
        
        return results