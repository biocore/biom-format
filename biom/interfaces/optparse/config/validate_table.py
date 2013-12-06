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
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.optparse.output_handler import print_list_of_strings
from biom.commands.table_validator import CommandConstructor
from biom.interfaces.optparse.input_handler import load_json_document

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

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
    OptparseOption(Parameter=cmd_in_lookup('table_json'),
                   Type='existing_filepath',
                   Handler=load_json_document, ShortName='i',
                   Name='input-fp',
                   Help='the input filepath to validate against the BIOM '
                   'format specification'),

    OptparseOption(Parameter=cmd_in_lookup('format_version'), ShortName='f'),

    OptparseOption(Parameter=cmd_in_lookup('detailed_report'), Type=None,
                   Action='store_true')
]

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('report_lines'),
                   Handler=print_list_of_strings),
    OptparseResult(Parameter=cmd_out_lookup('valid_table'),
                   Handler=report_table_validity)
]
