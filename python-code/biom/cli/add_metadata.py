#!/usr/bin/env python

from pyqi.interface.cli import CLOption, UsageExample, ParameterConversion
from biom.command import CommandConstructor

# How you can use the command from the command line
usage_examples = [
    ]

# Parameter conversions tell the interface how to describe command line 
# options
param_conversions = {
    		'biom-table':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='biom-table',
					CLType=<class 'biom.table.Table'>),
		'sample-mapping':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='sample-mapping',
					CLType=<class 'biom.parse.MetadataMap'>),
		'observation-mapping':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='observation-mapping',
					CLType=<class 'biom.parse.MetadataMap'>),
		'sc-separated':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='sc-separated',
					CLType=<type 'list'>),
		'sc-pipe-separated':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='sc-pipe-separated',
					CLType=<type 'list'>),
		'int-fields':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='int-fields',
					CLType=<type 'list'>),
		'float-fields':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='float-fields',
					CLType=<type 'list'>),
		'observation-header':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='observation-header',
					CLType=<type 'list'>),
		'sample-header':ParameterConversion(ShortName="MUST BE DEFINED",
					LongName='sample-header',
					CLType=<type 'list'>),

    }

# The output map associated keys in the results returned from Command.run
# without output handlers
output_map = {
    }

# In case there are interface specific bits such as output files
additional_options = [
    ]

