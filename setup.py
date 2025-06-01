#!/usr/bin/env python

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
__version__ = "2.1.16-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


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
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: 3.11
    Programming Language :: Python :: 3.12
    Programming Language :: Python :: 3.13
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: OS Independent
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
    Operating System :: Microsoft :: Windows
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

install_requires = [
    "click",
    "numpy >= 1.9.2",
    "scipy >= 1.8.0",
    'pandas >= 0.20.0',
    "h5py",
]

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
      tests_require=[
          'pytest>=6.2.4',
          'pytest-cov',
          'flake8',
      ],
      include_package_data=True,
      ext_modules=extensions,
      include_dirs=[np.get_include()],
      install_requires=install_requires,
      extras_require={'hdf5': ["h5py >= 2.2.0"],
                      'anndata': ["anndata"],
                      },
      classifiers=classifiers,
      entry_points='''
          [console_scripts]
          biom=biom.cli:cli
      ''')
