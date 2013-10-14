#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import sys
from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseResult)
from pyqi.core.command import make_parameter_collection_lookup_f
from pyqi.core.interfaces.optparse.output_handler import print_list_of_strings
from biom.commands.table_validator import CommandConstructor
from biom.interfaces.optparse.input_handler import load_json_document

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

param_lookup = make_parameter_collection_lookup_f(CommandConstructor)

def report_table_validity(result_key, data, option_value=None):
    if data:
        print "The input file is a valid BIOM-formatted file."
        sys.exit(0)
    else:
        print "The input file is not a valid BIOM-formatted file."
        sys.exit(1)

usage_examples = [
    OptparseUsageExample(ShortDesc="Validating a BIOM file",
                         LongDesc="Validate the contents of table.biom for "
                                  "adherence to the BIOM format specification",
                         Ex="%prog -i table.biom")
]

inputs = [
    OptparseOption(Parameter=param_lookup('table_json'),
                   InputType='existing_filepath',
                   InputHandler=load_json_document, ShortName='i',
                   Name='input-fp',
                   Help='the input filepath to validate against the BIOM '
                   'format specification'),

    OptparseOption(Parameter=param_lookup('format_version'), ShortName='f'),

    OptparseOption(Parameter=param_lookup('detailed_report'), InputType=None,
                   InputAction='store_true')
]

outputs = [
    OptparseResult(ResultKey='report_lines',
                   OutputHandler=print_list_of_strings),
    OptparseResult(ResultKey='valid_table',
                   OutputHandler=report_table_validity)
]
