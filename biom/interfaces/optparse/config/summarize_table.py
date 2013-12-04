#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout", "Daniel McDonald"]
__license__ = "BSD"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.optparse.output_handler import write_list_of_strings
from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseOption, OptparseResult)
from biom.commands.table_summarizer import CommandConstructor
from biom.interfaces.optparse.input_handler import (
        load_biom_table_with_file_contents)

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Basic script usage",
                         LongDesc="Write a summary of table.biom to table_summary.txt",
                         Ex="%prog -i table.biom -o table_summary.txt")
]

inputs = [
    OptparseOption(Parameter=cmd_in_lookup('table'),
                   Type="existing_filepath",
                   Handler=load_biom_table_with_file_contents,
                   ShortName='i',
                   Name='input-fp'),
    OptparseOption(Parameter=cmd_in_lookup('qualitative'),
                   Type=None,
                   Action="store_true"),
    OptparseOption(Parameter=cmd_in_lookup('suppress_md5'),
                   Type=None,
                   Action="store_true"),
    OptparseOption(Parameter=None,
                   Type='new_filepath',
                   ShortName='o',
                   Name='output-fp',
                   Required=True,
                   Help='the output filepath')
]

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('biom_summary'),
                   Handler=write_list_of_strings,
                   InputName='output-fp')
]
