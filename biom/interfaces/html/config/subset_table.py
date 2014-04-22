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
__credits__ = ["Evan Bolyen", "Jai Ram Rideout", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"

from pyqi.core.interfaces.html import (HTMLInputOption, HTMLDownload)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from pyqi.core.interfaces.html.input_handler import (load_file_contents,
                                                     load_file_lines)
from pyqi.core.interfaces.html.output_handler import newline_list_of_strings
from biom.commands.table_subsetter import CommandConstructor

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('table_str'),
                    Type='upload_file',
                    Handler=load_file_contents,
                    Name='input-fp',
                    Help='the input BIOM table file to subset'),

    HTMLInputOption(Parameter=cmd_in_lookup('axis')),

    HTMLInputOption(Parameter=cmd_in_lookup('ids'),
                    Type='upload_file', Handler=load_file_lines,
                    Help='a file containing a single column of '
                    'IDs to retain (either sample IDs or observation IDs, '
                    'depending on the axis)'),

    HTMLInputOption(Parameter=None,
                    Name='download-file',
                    Required=True,
                    Help='the download file')
]

outputs = [
    HTMLDownload(Parameter=cmd_out_lookup('subset_generator'),
                 Handler=newline_list_of_strings,
                 FilenameLookup='download-file',
                 FileExtension='.biom')
]
