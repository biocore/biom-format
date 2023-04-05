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
from biom.cli.util import write_biom_table
from biom.parse import MetadataMap


@cli.command(name='add-metadata')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('-o', '--output-fp', required=True,
              type=click.Path(exists=False, dir_okay=False),
              help='The output BIOM table')
@click.option('-m', '--sample-metadata-fp', required=False,
              type=click.Path(exists=True, dir_okay=False),
              help='The sample metadata mapping file (will add sample '
                   'metadata to the input BIOM table, if provided).')
@click.option('--observation-metadata-fp', required=False,
              type=click.Path(exists=True, dir_okay=False),
              help='The observation metadata mapping file (will add '
                   'observation metadata to the input BIOM table, if '
                   'provided).')
@click.option('--sc-separated', required=False, type=click.STRING,
              help='Comma-separated list of the metadata fields to split '
                   'on semicolons. This is useful for hierarchical data such '
                   'as taxonomy or functional categories.')
@click.option('--sc-pipe-separated', required=False, type=click.STRING,
              help='Comma-separated list of the metadata fields to split '
                   'on semicolons and pipes ("|"). This is useful for '
                   'hierarchical data such as functional categories with '
                   'one-to-many mappings (e.g. x;y;z|x;y;w)).')
@click.option('--int-fields', required=False, type=click.STRING,
              help='Comma-separated list of the metadata fields to cast '
                   'to integers. This is useful for integer data such as '
                   '"DaysSinceStart".')
@click.option('--float-fields', required=False, type=click.STRING,
              help='Comma-separated list of the metadata fields to cast '
                   'to floating point numbers. This is useful for real number '
                   'data such as "pH".')
@click.option('--sample-header', required=False, type=click.STRING,
              help='Comma-separated list of the sample metadata field '
                   'names. This is useful if a header line is not provided '
                   'with the metadata, if you want to rename the fields, or '
                   'if you want to include only the first n fields where n is '
                   'the number of entries provided here.')
@click.option('--observation-header', required=False, type=click.STRING,
              help='Comma-separated list of the observation metadata '
                   'field names. This is useful if a header line is not '
                   'provided with the metadata, if you want to rename the '
                   'fields, or if you want to include only the first n fields '
                   'where n is the number of entries provided here.')
@click.option('--output-as-json', default=False, is_flag=True,
              help='Write the output file in JSON format.')
def add_metadata(input_fp, output_fp, sample_metadata_fp,
                 observation_metadata_fp, sc_separated, sc_pipe_separated,
                 int_fields, float_fields, sample_header, observation_header,
                 output_as_json):
    """Add metadata to a BIOM table.

    Add sample and/or observation metadata to BIOM-formatted files. See
    examples here: http://biom-format.org/documentation/adding_metadata.html

    Example usage:

    Add sample metadata to a BIOM table:

    $ biom add-metadata -i otu_table.biom -o table_with_sample_metadata.biom
      -m sample_metadata.txt
    """
    table = load_table(input_fp)
    if sample_metadata_fp is not None:
        sample_metadata_f = open(sample_metadata_fp)
    else:
        sample_metadata_f = None
    if observation_metadata_fp is not None:
        observation_metadata_f = open(observation_metadata_fp)
    else:
        observation_metadata_f = None
    if sc_separated is not None:
        sc_separated = sc_separated.split(',')
    if sc_pipe_separated is not None:
        sc_pipe_separated = sc_pipe_separated.split(',')
    if int_fields is not None:
        int_fields = int_fields.split(',')
    if float_fields is not None:
        float_fields = float_fields.split(',')
    if sample_header is not None:
        sample_header = sample_header.split(',')
    if observation_header is not None:
        observation_header = observation_header.split(',')

    result = _add_metadata(table, sample_metadata_f, observation_metadata_f,
                           sc_separated, sc_pipe_separated, int_fields,
                           float_fields, sample_header, observation_header)

    if output_as_json:
        fmt = 'json'
    else:
        fmt = 'hdf5'

    write_biom_table(result, fmt, output_fp)


def _split_on_semicolons(x):
    return [e.strip() for e in x.split(';')]


def _split_on_semicolons_and_pipes(x):
    return [[e.strip() for e in y.split(';')] for y in x.split('|')]


def _int(x):
    try:
        return int(x)
    except ValueError:
        return x


def _float(x):
    try:
        return float(x)
    except ValueError:
        return x


def _add_metadata(table, sample_metadata=None, observation_metadata=None,
                  sc_separated=None, sc_pipe_separated=None, int_fields=None,
                  float_fields=None, sample_header=None,
                  observation_header=None):

    if sample_metadata is None and observation_metadata is None:
        raise ValueError('Must specify sample_metadata and/or '
                         'observation_metadata.')

    # define metadata processing functions, if any
    process_fns = {}
    if sc_separated is not None:
        process_fns.update(dict.fromkeys(sc_separated,
                                         _split_on_semicolons))

    if sc_pipe_separated is not None:
        process_fns.update(dict.fromkeys(sc_pipe_separated,
                           _split_on_semicolons_and_pipes))

    if int_fields is not None:
        process_fns.update(dict.fromkeys(int_fields, _int))

    if float_fields is not None:
        process_fns.update(dict.fromkeys(float_fields, _float))

    # parse mapping files
    if sample_metadata is not None:
        sample_metadata = MetadataMap.from_file(sample_metadata,
                                                process_fns=process_fns,
                                                header=sample_header)

    if observation_metadata is not None:
        observation_metadata = MetadataMap.from_file(
            observation_metadata,
            process_fns=process_fns,
            header=observation_header)

    # NAUGHTY: this is modifying the input table IN PLACE!!! And then
    # RETURNING IT! MetadataAdder is angry!

    # add metadata as necessary
    if sample_metadata:
        table.add_metadata(sample_metadata, axis='sample')

    if observation_metadata:
        table.add_metadata(observation_metadata, axis='observation')

    return table
