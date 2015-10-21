# ----------------------------------------------------------------------------
# Copyright (c) 2011-2015, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

from importlib import import_module

import click
import biom


@click.group()
@click.version_option(version=biom.__version__)
def cli():
    pass


import_module('biom.cli.table_summarizer')
import_module('biom.cli.metadata_adder')
import_module('biom.cli.table_converter')
import_module('biom.cli.installation_informer')
import_module('biom.cli.table_subsetter')
import_module('biom.cli.table_normalizer')
import_module('biom.cli.table_head')
import_module('biom.cli.table_validator')
import_module('biom.cli.uc_processor')
