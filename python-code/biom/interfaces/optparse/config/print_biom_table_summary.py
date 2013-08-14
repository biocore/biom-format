#!/usr/bin/env python

from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseOption, OptparseResult)
from pyqi.core.interfaces.optparse.output_handler import write_string
from biom.commands.summarize_biom_table import CommandConstructor
from biom.interfaces.optparse.input_handler import load_biom_table

usage_examples = [
    OptparseUsageExample(ShortDesc="A short single sentence description of the example",
                         LongDesc="A longer, more detailed description",
                         Ex="%prog --foo --bar some_file")
]

inputs = [
    OptparseOption(Parameter=CommandConstructor.Parameters['table'],
                   InputType="existing_filepath",
                   InputHandler=load_biom_table,
                   ShortName='i',
                   Name='input_fp',
                   Required=True,
                   Help='the input BIOM table'),
    OptparseOption(Parameter=CommandConstructor.Parameters['qualitative'],
                   InputType=None,
                   InputAction="store_true",
                   InputHandler=None,
                   ShortName=None,
                   Name='num_observations',
                   Required=False,
                   Help='Present counts as number of unique observation ids per sample, rather than counts of observations per sample.'),
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
                   OutputHandler=write_string,
                   OptionName='output_fp')
]
