#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The pyqi project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Doug Wendel",
               "Jai Ram Rideout"]
__license__ = "BSD"
__version__ = "0.1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

import os
import glob

# from http://stackoverflow.com/questions/1057431/loading-all-modules-in-a-folder-in-python
__all__ = sorted([os.path.basename(f)[:-3]
                  for f in glob.glob(os.path.dirname(__file__)+"/*.py")
                  if not os.path.basename(f).startswith('__init__')])
__all_lookup__ = set(__all__)