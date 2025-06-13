#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2020, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup
from setuptools.extension import Extension

long_description = """BIOM: Biological Observation Matrix
http://www.biom-format.org

The Biological Observation Matrix (BIOM) format or: how I learned to stop
worrying and love the ome-ome
Daniel McDonald, Jose C Clemente, Justin Kuczynski, Jai Ram Rideout,
Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle,
Folker Meyer, Rob Knight, J Gregory Caporaso
GigaScience 2012, 1:7.
"""

def get_extensions():
    import numpy as np
    from Cython.Build import cythonize

    extensions = [
        Extension("biom._filter",
                    ["biom/_filter.pyx"],
                    include_dirs=[np.get_include()]),
        Extension("biom._transform",
                    ["biom/_transform.pyx"],
                    include_dirs=[np.get_include()]),
        Extension("biom._subsample",
                    ["biom/_subsample.pyx"],
                    include_dirs=[np.get_include()]),
    ]
    return cythonize(extensions)


setup(ext_modules=get_extensions())