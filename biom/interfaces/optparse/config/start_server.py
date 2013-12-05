#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


__author__ = "Evan Bolyen"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Evan Bolyen"]
__license__ = "BSD"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"

from pyqi.core.interfaces.optparse import (OptparseOption,
                                           OptparseResult,
                                           OptparseUsageExample)
from pyqi.core.interfaces.optparse.input_handler import string_list_handler
from pyqi.core.interfaces.optparse.output_handler import print_string
from pyqi.core.command import make_command_in_collection_lookup_f, make_command_out_collection_lookup_f
from biom.commands.start_server import CommandConstructor

cmdin_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmdout_lookup = make_command_out_collection_lookup_f(CommandConstructor)



usage_examples = [
    OptparseUsageExample(ShortDesc="Start Biom's webserver",
                         LongDesc="Starts Biom's webserver on the specified --port",
                         Ex='%prog -p 8080 -m pyqi.interfaces.html.config')
]

inputs = [
    OptparseOption(Parameter=cmdin_lookup('port'),
                   ShortName='p',
                   Type=int)
]

outputs = [
    OptparseResult(Parameter=cmdout_lookup('result'),
                   Handler=print_string)
]
