cdef extern from "sparsemat_lib.h" namespace "sparsemat":
    cdef cppclass SparseMatFloat:
        SparseMatFloat()
        void insert(int, int, float)
        float get(int, int)

    cdef cppclass SparseMatInt:
        SparseMatInt()
        void insert(int, int, int)
        int get(int, int)

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
