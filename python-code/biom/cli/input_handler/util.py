#!/usr/bin/env python

"""Command line interface input handlers

All input handlers must conform to the following function definittion

function(option_value)
"""

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from biom.parse import parse_biom_table

def biom_table_handler(option_value):
    """Open and parse a BIOM table filepath."""
    return parse_biom_table(open(option_value, 'U'))
