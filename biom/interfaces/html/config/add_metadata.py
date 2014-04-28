#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from pyqi.core.interfaces.html import (HTMLInputOption, HTMLDownload)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.html.output_handler import newline_list_of_strings
from pyqi.core.interfaces.optparse.input_handler import string_list_handler

from biom.interfaces.html.input_handler import load_biom_table
from biom.commands.metadata_adder import CommandConstructor

__author__ = "Evan Bolyen"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = [
    "Evan Bolyen", "Jai Ram Rideout", "Greg Caporaso", "Morgan Langille",
    "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)


inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('table'),
                    Type='upload_file',
                    Handler=load_biom_table,
                    Name='input-fp'),

    HTMLInputOption(Parameter=cmd_in_lookup('sample_metadata'),
                    Type='upload_file',
                    Name='sample-metadata-fp'),

    HTMLInputOption(Parameter=cmd_in_lookup('observation_metadata'),
                    Type='upload_file',
                    Name='observation-metadata-fp'),

    HTMLInputOption(Parameter=cmd_in_lookup('sc_separated'),
                    Handler=string_list_handler,
                    Help='comma-separated list of the metadata fields to '
                    'split on semicolons. This is useful for hierarchical '
                    'data such as taxonomy or functional categories'),

    HTMLInputOption(Parameter=cmd_in_lookup('sc_pipe_separated'),
                    Handler=string_list_handler,
                    Help='comma-separated list of the metadata fields to split'
                    ' on semicolons and pipes ("|"). This is useful for '
                    'hierarchical data such as functional categories with '
                    'one-to-many mappings (e.g. x;y;z|x;y;w)'),

    HTMLInputOption(Parameter=cmd_in_lookup('int_fields'),
                    Handler=string_list_handler,
                    Help='comma-separated list of the metadata fields to cast '
                    'to integers. This is useful for integer data such as '
                    '"DaysSinceStart"'),

    HTMLInputOption(Parameter=cmd_in_lookup('float_fields'),
                    Handler=string_list_handler,
                    Help='comma-separated list of the metadata fields to cast '
                    'to floating point numbers. This is useful for real number'
                    ' data such as "pH"'),

    HTMLInputOption(Parameter=cmd_in_lookup('sample_header'),
                    Handler=string_list_handler,
                    Help='comma-separated list of the sample metadata field '
                    'names. This is useful if a header line is not provided '
                    'with the metadata, if you want to rename the fields, or '
                    'if you want to include only the first n fields where n is'
                    ' the number of entries provided here'),

    HTMLInputOption(Parameter=cmd_in_lookup('observation_header'),
                    Handler=string_list_handler,
                    Help='comma-separated list of the observation metadata '
                    'field names. This is useful if a header line is not '
                    'provided with the metadata, if you want to rename the '
                    'fields, or if you want to include only the first n fields'
                    ' where n is the number of entries provided here'),
    HTMLInputOption(Parameter=None,
                    Name='download-file',
                    Required=True,
                    Help='the download file')
]

outputs = [
    HTMLDownload(Parameter=cmd_out_lookup('table'),
                 Handler=newline_list_of_strings,
                 FilenameLookup='download-file',
                 FileExtension='.biom')
]
