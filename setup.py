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

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: OS Independent
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

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
    packages=['biom',
              'biom/backends',
              'biom/commands',
              'biom/interfaces',
              'biom/interfaces/optparse',
              'biom/interfaces/optparse/config'
              ],
    scripts=glob('scripts/*'),
    install_requires=["numpy >= 1.3.0",
                      "pyqi == 0.3.1"],
    extras_require={'scipy_sparse':["scipy >= 0.9.0"],
                    'test':["nose >= 0.10.1",
                            "tox >= 1.6.1"],
                    'validator':['dateutil >= 2.1']
                   },
    classifiers=classifiers
)
