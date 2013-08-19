#!/usr/bin/env python

from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseResult)
from pyqi.core.command import make_parameter_collection_lookup_f
from pyqi.core.interfaces.optparse.input_handler import (file_reading_handler,
                                                         string_list_handler)
from biom.commands.metadata_adder import CommandConstructor
from biom.interfaces.optparse.input_handler import load_biom_table
from biom.interfaces.optparse.output_handler import write_biom_table

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso", "Morgan Langille"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.2-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

param_lookup = make_parameter_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Adding sample metadata",
                         LongDesc="Add sample metadata to a BIOM table",
                         Ex="%prog -i otu_table.biom -o table_with_sample_metadata.biom -m sample_metadata.txt")
]

inputs = [
    OptparseOption(Parameter=param_lookup('table'),
                   InputType='existing_filepath',
                   InputHandler=load_biom_table, ShortName='i',
                   Name='input-fp'),

    OptparseOption(Parameter=param_lookup('sample_metadata'),
                   InputType='existing_filepath',
                   InputHandler=file_reading_handler, ShortName='m',
                   Name='sample-metadata-fp'),

    OptparseOption(Parameter=param_lookup('observation_metadata'),
                   InputType='existing_filepath',
                   InputHandler=file_reading_handler,
                   Name='observation-metadata-fp'),

    OptparseOption(Parameter=param_lookup('sc_separated'),
                   InputHandler=string_list_handler,
                   Help='comma-separated list of the metadata fields to split '
                   'on semicolons. This is useful for hierarchical data such '
                   'as taxonomy or functional categories'),

    OptparseOption(Parameter=param_lookup('sc_pipe_separated'),
                   InputHandler=string_list_handler,
                   Help='comma-separated list of the metadata fields to split '
                   'on semicolons and pipes ("|"). This is useful for '
                   'hierarchical data such as functional categories with '
                   'one-to-many mappings (e.g. x;y;z|x;y;w)'),

    OptparseOption(Parameter=param_lookup('int_fields'),
                   InputHandler=string_list_handler,
                   Help='comma-separated list of the metadata fields to cast '
                   'to integers. This is useful for integer data such as '
                   '"DaysSinceStart"'),

    OptparseOption(Parameter=param_lookup('float_fields'),
                   InputHandler=string_list_handler,
                   Help='comma-separated list of the metadata fields to cast '
                   'to floating point numbers. This is useful for real number '
                   'data such as "pH"'),

    OptparseOption(Parameter=param_lookup('sample_header'),
                   InputHandler=string_list_handler,
                   Help='comma-separated list of the sample metadata field '
                   'names. This is useful if a header line is not provided '
                   'with the metadata, if you want to rename the fields, or '
                   'if you want to include only the first n fields where n is '
                   'the number of entries provided here'),

    OptparseOption(Parameter=param_lookup('observation_header'),
                   InputHandler=string_list_handler,
                   Help='comma-separated list of the observation metadata '
                   'field names. This is useful if a header line is not '
                   'provided with the metadata, if you want to rename the '
                   'fields, or if you want to include only the first n fields '
                   'where n is the number of entries provided here'),

    OptparseOption(Parameter=None, InputType='new_filepath', ShortName='o',
                   Name='output-fp', Required=True,
                   Help='the output BIOM table')
]

outputs = [
    OptparseResult(ResultKey='table', OutputHandler=write_biom_table,
                   OptionName='output-fp')
]
