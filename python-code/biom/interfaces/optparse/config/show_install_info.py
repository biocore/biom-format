#!/usr/bin/env python

from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseResult)
from pyqi.core.command import make_parameter_collection_lookup_f
from pyqi.core.interfaces.optparse.output_handler import print_list_of_strings
from biom.commands.installation_informer import CommandConstructor

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.2.0"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

param_lookup = make_parameter_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Displaying installation info",
                         LongDesc="Display biom-format installation "
                                  "information",
                         Ex="%prog")
]

inputs = []

outputs = [
    OptparseResult(ResultKey='install_info_lines',
                   OutputHandler=print_list_of_strings)
]
