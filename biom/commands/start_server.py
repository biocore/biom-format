#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division

__credits__ = ["Evan Bolyen"]

from pyqi.core.command import (Command, CommandIn, CommandOut, ParameterCollection)
from pyqi.core.interfaces.html import start_server

class ServeHTMLInterface(Command):
    BriefDescription = "Start Biom's webserver"
    LongDescription = "Start Biom webserver on the specified port"
    CommandIns = ParameterCollection([
        CommandIn(Name='port', DataType=int,
                  Description='The port to run the server on', Required=False,
                  Default=8080)
        ])

    CommandOuts = ParameterCollection([
          CommandOut(Name='result',DataType=str, 
                    Description='Signals the termination of the HTMLInterface server')
          ])

    def run(self, **kwargs):
        """Start the HTMLInterface server with the port and interface_module"""
        fin = start_server(kwargs['port'], "biom.interfaces.html.config")

        return {'result': fin}

CommandConstructor = ServeHTMLInterface
