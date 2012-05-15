from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(ext_modules=[Extension(
                   "_sparsemat",                 # name of extension
                   ["sparsemat.pyx", "sparsemat_lib.cpp"], #  our Cython source
                   language="c++",  # causes Cython to create C++ source
                   library_dirs=['/Users/mcdonald/lib'],
                   include_dirs=['/Users/mcdonald/include/python2.7'])],
      cmdclass={'build_ext': build_ext})
