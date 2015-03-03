#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.html.output_handler import newline_list_of_strings
from pyqi.core.interfaces.html import (HTMLInputOption, HTMLDownload)
from biom.commands.table_summarizer import CommandConstructor
from biom.interfaces.html.input_handler import (
    load_biom_table_with_file_contents
    )

__author__ = "Evan Bolyen"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = [
    "Evan Bolyen",
    "Greg Caporaso",
    "Jai Ram Rideout",
    "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('table'),
                    Type="upload_file",
                    Handler=load_biom_table_with_file_contents,
                    Name='input-fp'),
    HTMLInputOption(Parameter=cmd_in_lookup('qualitative'),
                    Type=bool),
    HTMLInputOption(Parameter=None,
                    Name='download-file',
                    Required=True,
                    Help='the download file')
]

outputs = [
    HTMLDownload(Parameter=cmd_out_lookup('biom_summary'),
                 Handler=newline_list_of_strings,
                 FilenameLookup='download-file',
                 FileExtension='.biom.summary.txt')
]
