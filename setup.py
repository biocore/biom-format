#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2020, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys

from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.test import test as TestCommand
import numpy as np
from Cython.Build import cythonize

# Hack to prevent stupid "TypeError: 'NoneType' object is not callable" error
# in multiprocessing/util.py _exit_function when running `python
# setup.py test` (see
# http://www.eby-sarna.com/pipermail/peak/2010-May/003357.html),
# borrowed from https://github.com/getsentry/sentry/blob/master/setup.py
for m in ('multiprocessing', 'logging'):
    try:
        __import__(m)
    except ImportError:
        pass

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2020, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente",
               "Jai Ram Rideout", "Jorge Cañardo Alastuey", "Michael Hall"]
__license__ = "BSD"
__version__ = "2.1.10"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


# derived from https://docs.pytest.org/en/3.8.0/goodpractices.html
class PyTest(TestCommand):
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        try:
            import numpy
            try:
                # NumPy 1.14 changed repr output breaking our doctests,
                # request the legacy 1.13 style
                numpy.set_printoptions(legacy="1.13")
            except TypeError:
                # Old Numpy, output should be fine as it is :)
                # TypeError: set_printoptions() got an unexpected
                # keyword argument 'legacy'
                pass
        except ImportError:
            numpy = None

        import shlex

        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)


long_description = """BIOM: Biological Observation Matrix
http://www.biom-format.org

The Biological Observation Matrix (BIOM) format or: how I learned to stop
worrying and love the ome-ome
Daniel McDonald, Jose C Clemente, Justin Kuczynski, Jai Ram Rideout,
Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle,
Folker Meyer, Rob Knight, J Gregory Caporaso
GigaScience 2012, 1:7.
"""

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: 3.7
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: OS Independent
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

# Dealing with Cython
ext = '.pyx'
extensions = [Extension("biom._filter",
                        ["biom/_filter" + ext],
                        include_dirs=[np.get_include()]),
              Extension("biom._transform",
                        ["biom/_transform" + ext],
                        include_dirs=[np.get_include()]),
              Extension("biom._subsample",
                        ["biom/_subsample" + ext],
                        include_dirs=[np.get_include()])]
extensions = cythonize(extensions)

install_requires = ["click", "numpy >= 1.9.2", "future >= 0.16.0",
                    "scipy >= 1.3.1", 'pandas >= 0.20.0',
                    "six >= 1.10.0", "cython >= 0.29", "h5py",
                    "cython"]

if sys.version_info[0] < 3:
    raise SystemExit("Python 2.7 is no longer supported")


setup(name='biom-format',
      version=__version__,
      description='Biological Observation Matrix (BIOM) format',
      long_description=long_description,
      license=__license__,
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://www.biom-format.org',
      packages=find_packages(),
      tests_require=['pytest < 5.3.4',
                     'pytest-cov',
                     'flake8',
                     'nose'],
      include_package_data=True,
      ext_modules=extensions,
      include_dirs=[np.get_include()],
      install_requires=install_requires,
      extras_require={'hdf5': ["h5py >= 2.2.0"],
                      'anndata': ["anndata"],
                      },
      classifiers=classifiers,
      cmdclass={"pytest": PyTest},
      entry_points='''
          [console_scripts]
          biom=biom.cli:cli
      ''')
