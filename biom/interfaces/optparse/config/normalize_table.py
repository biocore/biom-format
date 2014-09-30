#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseUsageExample,
                                           OptparseResult)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from biom.interfaces.optparse.output_handler import write_biom_table
from biom.commands.table_normalizer import CommandConstructor

__author__ = "Michael Shaffer"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Michael Shaffer"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Michael Shaffer"
__email__ = "michael.shaffer@ucdenver.edu"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)


usage_examples = [
    OptparseUsageExample(ShortDesc="Normalizing a BIOM table to relative"
                                   "abundnace",
                         LongDesc="Take a BIOM table and replace all values "
                                  "with their relative abundance in relation "
                                  "to the sample",
                         Ex="%prog -i table.biom -r -o "
                            "normalized_table.biom"),
    OptparseUsageExample(ShortDesc="Converting a BIOM table to a "
                                   "presence/absence table",
                         LongDesc="Take a BIOM table and convert the values "
                                  "to 0's and 1's based on presensce or "
                                  "absence of observations",
                         Ex="%prog -i table.biom -p -o converted_table.biom")
]

inputs = [
    # table input
    OptparseOption(Parameter=cmd_in_lookup('biom_table'),
                   Type='existing_filepath',
                   Handler=None, ShortName='i',
                   Name='input-fp', Required=True,
                   Help='the input BIOM table filepath to subset'),

    # normalization to relative_abundance
    OptparseOption(Parameter=cmd_in_lookup('relative_abund'),
                   ShortName='r',
                   Action='store_true',
                   Help='convert table to relative abundance'),

    # conversion to presensce/absence
    OptparseOption(Parameter=cmd_in_lookup('presence_absence'),
                   ShortName='p',
                   Action='store_true',
                   Help='convert table to presence/absence'),

    # if relative abundance then normalize by sample or by observation
    OptparseOption(Parameter=cmd_in_lookup('axis'),
                   ShortName='a',
                   Default='sample'),

    OptparseOption(Parameter=None,
                   Type='new_filepath',
                   ShortName='o',
                   Name='output-fp',
                   Required=True,
                   Help='the output filepath')
]

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('table'),
                   Handler=write_biom_table,
                   InputName='output-fp')
]
