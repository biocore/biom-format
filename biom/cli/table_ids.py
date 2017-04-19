# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import click

from biom.cli import cli
from biom import load_table


@cli.command(name='table-ids')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('--observations', default=False, is_flag=True,
              help="Grab observation IDs")
def summarize_table(input_fp, observations):
    """Dump IDs in a table.

    Dump out the IDs found within a table:

    Example usage:

    Get the sample IDs within a table:

    $ biom table-ids -i table.biom

    Get the observation IDs within a table:

    $ biom table-ids -i table.biom --observations
    """
    tab = load_table(input_fp)
    for id_ in tab.ids(axis='observation' if observations else 'sample'):
        click.echo(id_)
