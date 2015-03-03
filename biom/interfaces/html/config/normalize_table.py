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

from biom.commands.table_converter import CommandConstructor

__author__ = "Michael Shaffer"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Michael Shaffer"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Michael Shaffer"
__email__ = "michael.shaffer@ucdenver.edu"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)


inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('biom_table'),
                    Type='upload_file',
                    Required=True,
                    Help='the input table filepath'),

    HTMLInputOption(Parameter=cmd_in_lookup('relative_abund'),
                    Type=bool),
    HTMLInputOption(Parameter=cmd_in_lookup('presence_absence'),
                    Type=bool),
    HTMLInputOption(Parameter=cmd_in_lookup('axis'),
                    Type='multiple_choice',
                    Choices=['observation', 'sample'],
                    Help='axis by which to normalize the table'),

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
