# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import click

from biom import load_table
from biom.cli import cli


@cli.command(name='export-metadata')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('-m', '--sample-metadata-fp', required=False,
              type=click.Path(exists=False, dir_okay=False),
              help='The sample metadata output file.')
@click.option('--observation-metadata-fp', required=False,
              type=click.Path(exists=False, dir_okay=False),
              help='The observation metadata output file.')
def export_metadata(input_fp, sample_metadata_fp, observation_metadata_fp):
    """Export metadata as TSV.

    Example usage:

    Export metadata as TSV:

    $ biom export-metadata -i otu_table.biom
      --sample-metadata-fp sample.tsv
      --observation-metadata-fp observation.tsv
    """
    table = load_table(input_fp)

    if sample_metadata_fp:
        _export_metadata(table, 'sample', input_fp, sample_metadata_fp)
    if observation_metadata_fp:
        _export_metadata(table, 'observation', input_fp,
                         observation_metadata_fp)


def _export_metadata(table, axis, input_fp, output_fp):
    try:
        metadata = table.metadata_to_dataframe(axis)
        metadata.to_csv(output_fp, sep='\t')
    except KeyError:
        click.echo('File {} does not contain {} metadata'.format(input_fp,
                                                                 axis))
