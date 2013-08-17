#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from sys import argv
from biom.util import old_to_new_biom_command

old_to_new_biom_command(argv)
