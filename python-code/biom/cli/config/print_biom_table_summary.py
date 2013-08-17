#!/usr/bin/env python

from pyqi.interface.cli import (CLOption,
                                UsageExample,
                                ParameterConversion,
                                OutputHandler)
from biom.command.print_biom_table_summary import CommandConstructor
from biom.cli.input_handler import biom_table_handler

def write_lines_to_filepath(result_key,data,option_value):
    if option_value is None:
        print '\n'.join(data)
    else:
        f = open(option_value,'w')
        f.write('\n'.join(data))
        f.close()

# How you can use the command from the command line
usage_examples = [
    UsageExample(ShortDesc="Print summary information on a biom file",
                 LongDesc="Print a summary on the observation counts per sample, including summary statistics, as well as metadata categories associated with samples and observations.",
                 Ex="%prog -i otu_table.biom -o table.biom")
    ]

# Parameter conversions tell the interface how to describe command line 
# options
param_conversions = {
        'table':ParameterConversion(ShortName="i",
                    LongName='input_fp',
                    CLType='existing_filepath',
                    InHandler=biom_table_handler,
                    ),
        'num-observations':ParameterConversion(ShortName=None,
                    LongName='num-observations',
                    CLType=None,
                    CLAction='store_true'),

    }

# The output map associated keys in the results returned from Command.run
# without output handlers
output_map = {'biom-summary':OutputHandler(OptionName='output_fp',
                                           Function=write_lines_to_filepath)
    }

# In case there are interface specific bits such as output files
additional_options = [
           CLOption(ShortName='o',
                    LongName='output_fp',
                    CLType='new_filepath',
                    Type=str,
                    Help='Path to write output summary',
                    Name='biom-summary',
                    Required=False,
                    Default=None,
                    DefaultDescription='written to stdout',
                    ResultName='biom-summary'),
           CLOption(LongName='suppress_md5',
                    CLType=None,
                    Type=bool,
                    Help=("Do not include the md5sum of the table in the "
                          "output (useful if you're concerned about runtime "
                          "of this script)"),
                    Name='suppress-md5',
                    Required=False,
                    Default=False),
    ]

