#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseResult)
from pyqi.core.command import make_command_out_collection_lookup_f
from pyqi.core.interfaces.optparse.output_handler import print_list_of_strings
from biom.commands.installation_informer import CommandConstructor

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Displaying installation info",
                         LongDesc="Display biom-format installation "
                                  "information",
                         Ex="%prog")
]

inputs = []

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('install_info_lines'),
                   Handler=print_list_of_strings)
]
