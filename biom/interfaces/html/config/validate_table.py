#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from pyqi.core.interfaces.html import (HTMLInputOption, HTMLPage)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from biom.commands.table_validator import CommandConstructor
from biom.interfaces.html.input_handler import load_json_document

__author__ = "Evan Bolyen"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Evan Bolyen", "Jai Ram Rideout", "Daniel McDonald",
               "Jorge Cañardo Alastuey"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)


def display_table_validity(result_key, data, option_value=None):
    if data is None:
        return "The input file is a valid BIOM-formatted file."
    else:
        to_join = ["The input file is not a valid BIOM-formatted file."]
        to_join += data
        return "<br/>".join(to_join)


inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('table'),
                    Type='upload_file',
                    Handler=load_json_document,
                    Name='input-fp',
                    Help='the input file to validate against the BIOM '
                    'format specification'),

    HTMLInputOption(Parameter=cmd_in_lookup('format_version')),

    HTMLInputOption(Parameter=cmd_in_lookup('detailed_report'), Type=bool),
]

outputs = [
    HTMLPage(Parameter=cmd_out_lookup('report_lines'),
             Handler=display_table_validity)
]
