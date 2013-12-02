#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from distutils.core import setup
from glob import glob
import numpy
import pyqi

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente",
               "Jai Ram Rideout"]
__license__ = "BSD"
__version__ = "1.2.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

long_description = """BIOM: Biological Observation Matrix
http://www.biom-format.org

The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome
Daniel McDonald, Jose C Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, J Gregory Caporaso
GigaScience 2012, 1:7.
"""

setup(name='biom-format',
    version=__version__,
    description='Biological Observation Matrix (BIOM) format',
    long_description=long_description,
    author=__maintainer__,
    author_email=__email__,
    maintainer=__maintainer__,
    maintainer_email=__email__,
    url='http://www.biom-format.org',
    packages=['biom',
              'biom/backends',
              'biom/commands',
              'biom/interfaces',
              'biom/interfaces/optparse',
              'biom/interfaces/optparse/config'
              ],
    scripts=glob('scripts/*')
)
