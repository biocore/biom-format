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
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.optparse.input_handler import (file_reading_handler,
                                                         string_list_handler)
from biom.commands.metadata_adder import CommandConstructor
from biom.interfaces.optparse.input_handler import load_biom_table
from biom.interfaces.optparse.output_handler import write_biom_table

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso", "Morgan Langille",
               "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Adding sample metadata",
                         LongDesc="Add sample metadata to a BIOM table",
                         Ex="%prog -i otu_table.biom -o table_with_sample_metadata.biom -m sample_metadata.txt")
]

inputs = [
    OptparseOption(Parameter=cmd_in_lookup('table'),
                   Type='existing_filepath',
                   Handler=load_biom_table, ShortName='i',
                   Name='input-fp'),

    OptparseOption(Parameter=cmd_in_lookup('sample_metadata'),
                   Type='existing_filepath',
                   Handler=file_reading_handler, ShortName='m',
                   Name='sample-metadata-fp'),

    OptparseOption(Parameter=cmd_in_lookup('observation_metadata'),
                   Type='existing_filepath',
                   Handler=file_reading_handler,
                   Name='observation-metadata-fp'),

    OptparseOption(Parameter=cmd_in_lookup('sc_separated'),
                   Handler=string_list_handler,
                   Help='comma-separated list of the metadata fields to split '
                   'on semicolons. This is useful for hierarchical data such '
                   'as taxonomy or functional categories'),

    OptparseOption(Parameter=cmd_in_lookup('sc_pipe_separated'),
                   Handler=string_list_handler,
                   Help='comma-separated list of the metadata fields to split '
                   'on semicolons and pipes ("|"). This is useful for '
                   'hierarchical data such as functional categories with '
                   'one-to-many mappings (e.g. x;y;z|x;y;w)'),

    OptparseOption(Parameter=cmd_in_lookup('int_fields'),
                   Handler=string_list_handler,
                   Help='comma-separated list of the metadata fields to cast '
                   'to integers. This is useful for integer data such as '
                   '"DaysSinceStart"'),

    OptparseOption(Parameter=cmd_in_lookup('float_fields'),
                   Handler=string_list_handler,
                   Help='comma-separated list of the metadata fields to cast '
                   'to floating point numbers. This is useful for real number '
                   'data such as "pH"'),

    OptparseOption(Parameter=cmd_in_lookup('sample_header'),
                   Handler=string_list_handler,
                   Help='comma-separated list of the sample metadata field '
                   'names. This is useful if a header line is not provided '
                   'with the metadata, if you want to rename the fields, or '
                   'if you want to include only the first n fields where n is '
                   'the number of entries provided here'),

    OptparseOption(Parameter=cmd_in_lookup('observation_header'),
                   Handler=string_list_handler,
                   Help='comma-separated list of the observation metadata '
                   'field names. This is useful if a header line is not '
                   'provided with the metadata, if you want to rename the '
                   'fields, or if you want to include only the first n fields '
                   'where n is the number of entries provided here'),

    OptparseOption(Parameter=None, Type='new_filepath', ShortName='o',
                   Name='output-fp', Required=True,
                   Help='the output BIOM table')
]

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('table'),
                   Handler=write_biom_table,
                   InputName='output-fp')
]
