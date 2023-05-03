#!/usr/bin/env python

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
from biom.cli.util import write_biom_table


@cli.command(name='normalize-table')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('-o', '--output-fp', default=None,
              type=click.Path(writable=True),
              help='An output file-path')
@click.option('-r', '--relative-abund', default=False, is_flag=True,
              help='convert table to relative abundance',
              required=False)
@click.option('-p', '--presence-absence', default=False, is_flag=True,
              help='convert table to presence/absence',
              required=False)
@click.option('-a', '--axis', default='sample',
              type=click.Choice(['sample', 'observation']),
              help='The axis to normalize over')
def normalize_table(input_fp, output_fp, relative_abund, presence_absence,
                    axis):
    """Normalize a BIOM table.

    Normalize the values of a BIOM table through various methods. Relative
    abundance will take the relative abundance of each observation in terms of
    samples or observations.  Presence absensece will convert observations to
    1's and 0's based on presence of the observation.

    Example usage:

    Normalizing a BIOM table to relative abundnace:

    $ biom normalize-table -i table.biom -r -o normalized_table.biom

    Converting a BIOM table to a presence/absence table:

    $ biom normalize-table -i table.biom -p -o converted_table.biom
    """
    table = load_table(input_fp)
    result = _normalize_table(table, relative_abund, presence_absence, axis)

    write_biom_table(result, 'hdf5', output_fp)


def _normalize_table(table, relative_abund=False, presence_absence=False,
                     axis='sample'):
    if relative_abund is False and presence_absence is False:
        raise ValueError("Must specifiy a normalization type")
    elif relative_abund is True and presence_absence is True:
        raise ValueError("Must specify only one normalization type")

    if relative_abund:
        table.norm(axis=axis)
    else:
        table.pa()

    return table
