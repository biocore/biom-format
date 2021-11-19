# ----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from importlib import import_module

import click
import biom


def _terribly_handle_brokenpipeerror():
    # based off http://stackoverflow.com/a/34299346
    import os
    import sys
    sys.stdout = os.fdopen(1, 'w')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version=biom.__version__)
@click.pass_context
def cli(ctx):
    ctx.call_on_close(_terribly_handle_brokenpipeerror)


import_module('biom.cli.table_summarizer')
import_module('biom.cli.metadata_adder')
import_module('biom.cli.metadata_exporter')
import_module('biom.cli.table_converter')
import_module('biom.cli.installation_informer')
import_module('biom.cli.table_subsetter')
import_module('biom.cli.table_normalizer')
import_module('biom.cli.table_head')
import_module('biom.cli.table_ids')
import_module('biom.cli.table_validator')
import_module('biom.cli.uc_processor')
