#!/usr/bin/env python

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from pyqi.core.interfaces.optparse import (OptparseUsageExample,
                                           OptparseOption, OptparseResult)
from pyqi.core.command import make_parameter_collection_lookup_f
from pyqi.core.interfaces.optparse.input_handler import file_reading_handler
from pyqi.core.interfaces.optparse.output_handler import write_string
from biom.interfaces.optparse.input_handler import load_metadata
from biom.commands.table_converter import CommandConstructor

param_lookup = make_parameter_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Converting from classic to BIOM format",
                         LongDesc="Convert the classic file table.txt to a "
                                  "sparse BIOM format OTU table",
                         Ex='%prog -i table.txt -o table.biom '
                            '--table-type "otu table"')
]

inputs = [
    OptparseOption(Parameter=param_lookup('table_file'),
                   InputType='existing_filepath',
                   InputHandler=file_reading_handler,
                   ShortName='i', Name='input-fp',
                   Help='the input table filepath, either in BIOM or classic '
                   'format'),
    OptparseOption(Parameter=param_lookup('matrix_type'), ShortName='t'),
    OptparseOption(Parameter=param_lookup('biom_to_classic_table'),
                   InputType=None, InputAction='store_true', ShortName='b'),
    OptparseOption(Parameter=param_lookup('sparse_biom_to_dense_biom'),
                   InputType=None,
                   InputAction='store_true'),
    OptparseOption(Parameter=param_lookup('dense_biom_to_sparse_biom'),
                   InputType=None,
                   InputAction='store_true'),
    OptparseOption(Parameter=param_lookup('sample_metadata'),
                   InputType='existing_filepath',
                   InputHandler=load_metadata,
                   ShortName='m',
                   Name='sample-metadata-fp'),
    OptparseOption(Parameter=param_lookup('observation_metadata'),
                   InputType='existing_filepath',
                   InputHandler=load_metadata, Name='observation-metadata-fp'),
    OptparseOption(Parameter=param_lookup('header_key')),
    OptparseOption(Parameter=param_lookup('output_metadata_id')),
    OptparseOption(Parameter=param_lookup('process_obs_metadata')),
    OptparseOption(Parameter=param_lookup('table_type')),
    OptparseOption(Parameter=None,
                   InputType='new_filepath',
                   ShortName='o',
                   Name='output-fp',
                   Required=True,
                   Help='the output filepath')
]

outputs = [
    OptparseResult(ResultKey='table_str', OutputHandler=write_string,
                   OptionName='output-fp')
]
