#!/usr/bin/env python

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

import os
import glob

# from http://stackoverflow.com/questions/1057431/loading-all-modules-in-a-folder-in-python
__all__ = [os.path.basename(f)[:-3] for f in glob.glob(os.path.dirname(__file__)+"/*.py")
                             if not os.path.basename(f).startswith('__init__')]
__all_lookup__ = set(__all__)
