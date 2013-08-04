#!/usr/bin/env python

from pyqi.interface.cli import CLOption, UsageExample, ParameterConversion
from pyqi.interface.input_handler.cli import file_reading_handler

from biom.command.add_metadata import CommandConstructor
from biom.cli.input_handler import biom_table_handler, string_list_handler

usage_examples = [
    UsageExample(ShortDesc="Adding sample metadata",
                 LongDesc="Add sample metadata to a BIOM table",
                 Ex="%prog -i otu_table.biom -o table-with-sample-metadata.biom --sample_mapping_fp sample_metadata.txt")
]

param_conversions = {
    'biom-table': ParameterConversion(ShortName='i', LongName='input_fp',
                                      CLType='existing_filepath',
                                      InHandler=biom_table_handler),
    'sample-mapping': ParameterConversion(ShortName='m',
                                          LongName='sample_mapping_fp',
                                          CLType='existing_filepath',
                                          InHandler=file_reading_handler),
    'observation-mapping': ParameterConversion(ShortName=None,
            LongName='observation_mapping_fp', CLType='existing_filepath',
            InHandler=file_reading_handler),
    'sc-separated': ParameterConversion(ShortName=None,
                                        LongName='sc_separated',
                                        CLType=str,
                                        InHandler=string_list_handler),
    'sc-pipe-separated': ParameterConversion(ShortName=None,
                                             LongName='sc_pipe_separated',
                                             CLType=str,
                                             InHandler=string_list_handler),
    'int-fields': ParameterConversion(ShortName=None,
                                      LongName='int_fields', CLType=str,
                                      InHandler=string_list_handler),
    'float-fields': ParameterConversion(ShortName=None,
                                        LongName='float_fields', CLType=str,
                                        InHandler=string_list_handler),
    'observation-header': ParameterConversion(ShortName=None,
                                              LongName='observation_header',
                                              CLType=str,
                                              InHandler=string_list_handler),
    'sample-header': ParameterConversion(ShortName=None,
                                         LongName='sample_header', CLType=str,
                                         InHandler=string_list_handler)
}

output_map = {
    'biom-table': OutputHandler(OutputName='output_fp',
                                Function=write_biom_table)
}

additional_options = [
    CLOption(Type=Table,
             Help='the output BIOM table',
             Name='biom-table',
             Required=True,
             LongName='output_fp',
             CLType='new_filepath',
             ShortName='o')
]
