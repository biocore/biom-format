from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy

# Compiling Pyrex modules to .c and .so
include_path = os.path.join(os.getcwd(), 'include')
# find arrayobject.h on every system an alternative would be to put
# arrayobject.h into pycogent/include, but why .. 
library_path = os.path.split(numpy.__file__)[0]

setup(ext_modules=[Extension(
                   "_sparsemat",                 # name of extension
                   ["sparsemat.pyx", "sparsemat_lib.cpp"], #  our Cython source
                   language="c++",  # causes Cython to create C++ source
                   library_dirs=library_path,
                   include_dirs=include_path)],
      cmdclass={'build_ext': build_ext})
