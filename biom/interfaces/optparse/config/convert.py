#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from pyqi.core.interfaces.optparse import (OptparseUsageExample,
                                           OptparseOption, OptparseResult)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.optparse.input_handler import file_reading_handler
from pyqi.core.interfaces.optparse.output_handler import write_string
from biom.interfaces.optparse.input_handler import load_metadata
from biom.commands.table_converter import CommandConstructor

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Converting from classic to BIOM format",
                         LongDesc="Convert the classic file table.txt to a "
                                  "sparse BIOM format OTU table",
                         Ex='%prog -i table.txt -o table.biom '
                            '--table-type "otu table"')
]

inputs = [
    OptparseOption(Parameter=cmd_in_lookup('table_file'),
                   Type='existing_filepath',
                   Handler=file_reading_handler,
                   ShortName='i', Name='input-fp',
                   Help='the input table filepath, either in BIOM or classic '
                   'format'),
    OptparseOption(Parameter=cmd_in_lookup('matrix_type'), ShortName='t'),
    OptparseOption(Parameter=cmd_in_lookup('biom_to_classic_table'),
                   Type=None, Action='store_true', ShortName='b'),
    OptparseOption(Parameter=cmd_in_lookup('sparse_biom_to_dense_biom'),
                   Type=None,
                   Action='store_true'),
    OptparseOption(Parameter=cmd_in_lookup('dense_biom_to_sparse_biom'),
                   Type=None,
                   Action='store_true'),
    OptparseOption(Parameter=cmd_in_lookup('sample_metadata'),
                   Type='existing_filepath',
                   Handler=load_metadata,
                   ShortName='m',
                   Name='sample-metadata-fp'),
    OptparseOption(Parameter=cmd_in_lookup('observation_metadata'),
                   Type='existing_filepath',
                   Handler=load_metadata, Name='observation-metadata-fp'),
    OptparseOption(Parameter=cmd_in_lookup('header_key')),
    OptparseOption(Parameter=cmd_in_lookup('output_metadata_id')),
    OptparseOption(Parameter=cmd_in_lookup('process_obs_metadata')),
    OptparseOption(Parameter=cmd_in_lookup('table_type')),
    OptparseOption(Parameter=None,
                   Type='new_filepath',
                   ShortName='o',
                   Name='output-fp',
                   Required=True,
                   Help='the output filepath')
]

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('table_str'), 
                   Handler=write_string,
                   InputName='output-fp')
]
