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
                   Name='input_fp', convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('sample_metadata'),
                   InputType='existing_filepath',
                   InputHandler=file_reading_handler, ShortName='m',
                   Name='sample_mapping_fp', convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('observation_metadata'),
                   InputType='existing_filepath',
                   InputHandler=file_reading_handler,
                   Name='observation_mapping_fp',
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('sc_separated'),
                   InputHandler=string_list_handler,
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('sc_pipe_separated'),
                   InputHandler=string_list_handler,
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('int_fields'),
                   InputHandler=string_list_handler,
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('float_fields'),
                   InputHandler=string_list_handler,
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('sample_header'),
                   InputHandler=string_list_handler,
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=param_lookup('observation_header'),
                   InputHandler=string_list_handler,
                   convert_to_dashed_name=False),

    OptparseOption(Parameter=None, InputType='new_filepath', ShortName='o',
                   Name='output_fp', Required=True,
                   Help='the output BIOM table', convert_to_dashed_name=False)
]

outputs = [
    OptparseResult(ResultKey='table', OutputHandler=write_biom_table,
                   OptionName='output_fp')
]
