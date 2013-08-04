#!/usr/bin/env python

"""Command line interface output handlers

All output handlers must conform to the following function definition

function(result_key, data, option_value=None)

result_key   - the corresponding key in the results dictionary
data         - the actual results
option_value - if the handler is tied to an output option, the value of that
               option is represented here
"""

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from pyqi.core.exception import IncompetentDeveloperError
import os
from biom.parse import generatedby

def write_biom_table(result_key, data, option_value=None):
    """Write a string to a file"""
    if option_value is None:
        raise IncompetentDeveloperError("Cannot write output without a "
                                        "filepath.")

    if os.path.exists(option_value):
        raise IOError("Output path '%s' already exists." % option_value)

    f = open(option_value, 'w')
    f.write(data.getBiomFormatJsonString(generatedby()))
    f.close()
