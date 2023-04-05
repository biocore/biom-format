# ----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


import click

from biom import load_table
from biom.cli import cli


@cli.command()
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('-o', '--output-fp', default=None,
              type=click.Path(writable=True),
              help='An output file-path', required=False)
@click.option('-n', '--n-obs', default=5, type=int,
              help="The number of observations to show",
              required=False)
@click.option('-m', '--n-samp', default=5, type=int,
              help="The number of samples to show",
              required=False)
def head(input_fp, output_fp, n_obs, n_samp):
    """Dump the first bit of a table.

    Example usage:

    Print out the upper left corner of a BIOM table to standard out:

    $ biom head -i table.biom

    """
    if n_obs == 0 or n_samp == 0:
        raise ValueError("Sample and observation identifiers can be obtained "
                         "from 'biom table-ids'")

    if n_obs < 0:
        raise ValueError("-n/--n-obs must be > 0")

    if n_samp < 0:
        raise ValueError("-m/--m-samp must be > 0")

    table = load_table(input_fp).head(n=n_obs, m=n_samp)

    if output_fp is None:
        click.echo(str(table))
    else:
        with open(output_fp, 'w') as fp:
            fp.write(str(table))
