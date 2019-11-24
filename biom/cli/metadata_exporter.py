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
def add_metadata(input_fp, sample_metadata_fp, observation_metadata_fp):
    """Export metadata as TSV.

    Example usage:

    Export metadata as TSV:

    $ biom export-metadata -i otu_table.biom
      --sample-metadata-fp sample.tsv
      --observation-metadata-fp observation.tsv
    """
    table = load_table(input_fp)
    sample_metadata = table.metadata_to_dataframe('sample')
    observation_metadata = table.metadata_to_dataframe('observation')

    if sample_metadata_fp:
        sample_metadata.to_csv(sample_metadata_fp, sep='\t')
    if observation_metadata_fp:
        observation_metadata.to_csv(observation_metadata_fp, sep='\t')
