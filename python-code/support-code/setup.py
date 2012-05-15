from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from os.path import join, split
from os import getcwd
try:
    import numpy
except ImportError:
    raise ImportError, "numpy cannot be found. Can't continue."

include_path = join(getcwd(), 'include')
library_path = split(numpy.__file__)[0]

setup(ext_modules=[Extension(
                   "_sparsemat",                 # name of extension
                   ["sparsemat.pyx", "sparsemat_lib.cpp"], #  our Cython source
                   language="c++",  # causes Cython to create C++ source
                   library_dirs=[library_path],
                   include_dirs=[include_path])],
      cmdclass={'build_ext': build_ext})
