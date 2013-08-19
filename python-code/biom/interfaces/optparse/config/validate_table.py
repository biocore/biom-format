#!/usr/bin/env python

import sys
from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseResult)
from pyqi.core.command import make_parameter_collection_lookup_f
from pyqi.core.interfaces.optparse.input_handler import file_reading_handler
from pyqi.core.interfaces.optparse.output_handler import print_list_of_strings
from biom.commands.table_validator import CommandConstructor

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.2-dev"
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
    OptparseOption(Parameter=param_lookup('table_file'),
                   InputType='existing_filepath',
                   InputHandler=file_reading_handler, ShortName='i',
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
