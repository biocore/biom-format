#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from pyqi.core.interfaces.html import HTMLPage
from pyqi.core.command import make_command_out_collection_lookup_f
from pyqi.core.interfaces.html.output_handler import html_list_of_strings
from biom.commands.installation_informer import CommandConstructor

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

cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)


inputs = []

outputs = [
    HTMLPage(Parameter=cmd_out_lookup('install_info_lines'),
             Handler=html_list_of_strings)
]
