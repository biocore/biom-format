#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

__author__ = "Evan Bolyen"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = [
    "Evan Bolyen",
    "Jai Ram Rideout",
    "Greg Caporaso",
    "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"

from pyqi.core.interfaces.html import (HTMLInputOption, HTMLDownload)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.html.output_handler import newline_list_of_strings

from biom.interfaces.html.input_handler import load_metadata
from biom.commands.table_converter import CommandConstructor

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)


inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('table_file'),
                    Type='upload_file',
                    Help='the input table filepath, either in BIOM or classic '
                    'format'),
    HTMLInputOption(Parameter=cmd_in_lookup('matrix_type'),
                    Type='multiple_choice',
                    Choices=['sparse', 'dense'],
                    Help='the type of BIOM file to create when a classic '
                    'table is supplied'),
    HTMLInputOption(Parameter=cmd_in_lookup('biom_to_classic_table'),
                    Type=bool),
    HTMLInputOption(Parameter=cmd_in_lookup('sparse_biom_to_dense_biom'),
                    Type=bool),
    HTMLInputOption(Parameter=cmd_in_lookup('dense_biom_to_sparse_biom'),
                    Type=bool),
    HTMLInputOption(Parameter=cmd_in_lookup('sample_metadata'),
                    Type='upload_file',
                    Handler=load_metadata),
    HTMLInputOption(Parameter=cmd_in_lookup('observation_metadata'),
                    Type='upload_file',
                    Handler=load_metadata),
    HTMLInputOption(Parameter=cmd_in_lookup('header_key')),
    HTMLInputOption(Parameter=cmd_in_lookup('output_metadata_id')),
    HTMLInputOption(Parameter=cmd_in_lookup('process_obs_metadata'),
                    Type='multiple_choice',
                    Choices=['taxonomy', 'naive', 'sc_separated'],
                    Help='Process metadata associated with observations when '
                    'converting from a classic table'),
    HTMLInputOption(Parameter=cmd_in_lookup('table_type'),
                    Type='multiple_choice',
                    Choices=[
                        'metabolite table',
                        'gene table',
                        'otu table',
                        'pathway table',
                        'function table',
                        'ortholog table',
                        'taxon table'],
                    Help='The BIOM table type to get converted into. Required '
                    'when converting a classic table file to a BIOM table '
                    'file.'),
    HTMLInputOption(Parameter=None,
                    Name='download-file',
                    Required=True,
                    Help='the download file')
]

outputs = [
    HTMLDownload(Parameter=cmd_out_lookup('table_str'),
                 Handler=newline_list_of_strings,
                 FilenameLookup='download-file',
                 FileExtension='.biom')
]
