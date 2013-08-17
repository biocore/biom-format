#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from pyqi.core.command import make_parameter_collection_lookup_f
from pyqi.core.interfaces.optparse.output_handler import write_list_of_strings
from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseOption, OptparseResult)
from biom.commands.table_summarizer import CommandConstructor
from biom.interfaces.optparse.input_handler import (
        load_biom_table_with_file_contents)

param_lookup = make_parameter_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Basic script usage",
                         LongDesc="Write a summary of table.biom to table_summary.txt",
                         Ex="%prog -i table.biom -o table_summary.txt")
]

inputs = [
    OptparseOption(Parameter=param_lookup('table'),
                   InputType="existing_filepath",
                   InputHandler=load_biom_table_with_file_contents,
                   ShortName='i',
                   Name='input_fp'),
    OptparseOption(Parameter=param_lookup('qualitative'),
                   InputType=None,
                   InputAction="store_true",
                   InputHandler=None,
                   ShortName=None,
                   Name='num_observations'),
    OptparseOption(Parameter=param_lookup('suppress_md5'),
                   InputType=None,
                   InputAction="store_true",
                   InputHandler=None,
                   ShortName=None,
                   Name='suppress_md5'),
    OptparseOption(Parameter=None,
                   InputType='new_filepath',
                   InputHandler=None,
                   ShortName='o',
                   Name='output_fp',
                   Required=True,
                   Help='output filepath')
]

outputs = [
    OptparseResult(ResultKey='biom-summary',
                   OutputHandler=write_list_of_strings,
                   OptionName='output_fp')
]
