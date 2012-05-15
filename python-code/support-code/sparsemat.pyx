# to int or to uint...
cdef extern from "sparsemat_lib.h" namespace "sparsemat":
    struct items_int:
        int *rows
        int *cols
        int *values
    
    struct items_float:
        int *rows
        int *cols
        float *values
        
cdef extern from "sparsemat_lib.h" namespace "sparsemat":
    cdef cppclass SparseMatFloat:
        SparseMatFloat()
        void insert(int, int, float)
        float get(int, int)
        void erase(int, int)
        int contains(int, int)
        int length()
        items_float keys()
        
    cdef cppclass SparseMatInt:
        SparseMatInt()
        void insert(int, int, int)
        int get(int, int)
        void erase(int, int)
        int contains(int, int)
        int length()
        
# should really use inheritence here but things were getting odd with *thisptr
# ...and i dont care right now

cdef class PySparseMatFloat:
    cdef SparseMatFloat *thisptr
    def __cinit__(self):
        self.thisptr = new SparseMatFloat()

    def __dealloc__(self):
        del self.thisptr

    def insert(self, row, col, value):
        self.thisptr.insert(<unsigned int>row, <unsigned int>col, value)
   
    def get(self, row, col):
        return self.thisptr.get(<unsigned int>row, <unsigned int>col)

    def erase(self, row, col):
        self.thisptr.erase(<unsigned int>row, <unsigned int>col)
        
    def contains(self, row, col):
        return self.thisptr.contains(<unsigned int>row, <unsigned int>col)
    
    def length(self):
        return self.thisptr.length()
        
    def keys(self):
        cdef items_float foo
        foo = self.thisptr.keys()
        
        results = []
        for i in range(self.length()):
            results.append((foo.rows[i], foo.cols[i]))
        
        return results
    
cdef class PySparseMatInt:
    cdef SparseMatInt *thisptr
    def __cinit__(self):
        self.thisptr = new SparseMatInt()

    def __dealloc__(self):
        del self.thisptr
   
    def insert(self, row, col, value):
        self.thisptr.insert(<unsigned int>row, <unsigned int>col, value)
   
    def get(self, row, col):
        return self.thisptr.get(<unsigned int>row, <unsigned int>col)

    def erase(self, row, col):
        self.thisptr.erase(<unsigned int>row, <unsigned int>col)

    def contains(self, row, col):
        return self.thisptr.contains(<unsigned int>row, <unsigned int>col)

    def length(self):
        return self.thisptr.length()
        
