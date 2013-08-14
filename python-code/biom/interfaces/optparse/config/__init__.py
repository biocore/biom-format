#!/usr/bin/env python

####
# This file is copied from pyqi/interfaces/optparse/config/__init__.py
# I think we want to work out a better solution for this as we'll need 
# this file a lot.
####

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The pyqi project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Doug Wendel",
               "Jai Ram Rideout"]
__license__ = "GPL"
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