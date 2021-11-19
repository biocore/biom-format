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


table_types = ["OTU table",
               "Pathway table",
               "Function table",
               "Ortholog table",
               "Gene table",
               "Metabolite table",
               "Taxon table",
               "Table"]

observation_metadata_types = {
    'sc_separated': lambda x: [e.strip() for e in x.split(';')],
    'naive': lambda x: x
}
observation_metadata_types['taxonomy'] = \
    observation_metadata_types['sc_separated']

observation_metadata_formatters = {
    'sc_separated': lambda x: '; '.join(x),
    'naive': lambda x: x
}


@cli.command(name='convert')
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
@click.option('--to-json', default=False, is_flag=True,
              help='Output as JSON-formatted table.')
@click.option('--to-hdf5', default=False, is_flag=True,
              help='Output as HDF5-formatted table.')
@click.option('--to-tsv', default=False, is_flag=True,
              help='Output as TSV-formatted (classic) table.')
@click.option('--collapsed-samples', default=False, is_flag=True,
              help='If --to_hdf5 is passed and the original table is a '
                   'BIOM table with collapsed samples, this will '
                   'update the sample metadata of the table to '
                   'the supported HDF5 collapsed format.')
@click.option('--collapsed-observations', default=False, is_flag=True,
              help='If --to_hdf5 is passed and the original table is a '
                   'BIOM table with collapsed observations, this will '
                   'update the observation metadata of the table '
                   'to the supported HDF5 collapsed format.')
@click.option('--header-key', required=False, type=click.STRING,
              help='The observation metadata to include from the input '
                   'BIOM table file when creating a tsv table file. '
                   'By default no observation metadata will be included.')
@click.option('--output-metadata-id', required=False, type=click.STRING,
              help='The name to be given to the observation metadata '
                   'column when creating a tsv table file if the column '
                   'should be renamed.')
@click.option('--table-type', required=False,
              type=click.Choice(table_types),
              help='The type of the table.')
@click.option('--process-obs-metadata', required=False,
              type=click.Choice(
                observation_metadata_types),
              help='Process metadata associated with observations when '
              'converting from a classic table.')
@click.option('--tsv-metadata-formatter', required=False,
              default='sc_separated',
              type=click.Choice(
                observation_metadata_formatters),
              help='Method for formatting the observation metadata.')
def convert(input_fp, output_fp, sample_metadata_fp, observation_metadata_fp,
            to_json, to_hdf5, to_tsv, collapsed_samples,
            collapsed_observations, header_key, output_metadata_id, table_type,
            process_obs_metadata, tsv_metadata_formatter):
    """Convert to/from the BIOM table format.

    Convert between BIOM table formats. See examples here:
    http://biom-format.org/documentation/biom_conversion.html

    Example usage:

    Convert a "classic" BIOM file (tab-separated text) to an HDF5 BIOM
    formatted OTU table:

    $ biom convert -i table.txt -o table.biom --to-hdf5
    """
    if sum([to_tsv, to_hdf5, to_json]) > 1:
        raise ValueError("--to-tsv, --to-json, and --to-hdf5 are mutually "
                         "exclusive. You can only pass one of these options.")

    table = load_table(input_fp)
    if sample_metadata_fp is not None:
        with open(sample_metadata_fp) as f:
            sample_metadata_f = MetadataMap.from_file(f)
    else:
        sample_metadata_f = None
    if observation_metadata_fp is not None:
        with open(observation_metadata_fp) as f:
            observation_metadata_f = MetadataMap.from_file(f)
    else:
        observation_metadata_f = None

    _convert(table, output_fp, sample_metadata_f, observation_metadata_f,
             to_json, to_hdf5, to_tsv, collapsed_samples,
             collapsed_observations, header_key, output_metadata_id,
             table_type, process_obs_metadata, tsv_metadata_formatter)


def _convert(table, output_filepath, sample_metadata=None,
             observation_metadata=None, to_json=False, to_hdf5=False,
             to_tsv=False, collapsed_samples=False,
             collapsed_observations=False, header_key=None,
             output_metadata_id=None, table_type=None,
             process_obs_metadata=None, tsv_metadata_formatter='sc_separated'):

    if sum([to_tsv, to_hdf5, to_json]) == 0:
        raise ValueError("Must specify an output format")
    elif sum([to_tsv, to_hdf5, to_json]) > 1:
        raise ValueError("Can only specify a single output format")

    if table_type is None:
        if table.type in [None, "None"]:
            table.type = "Table"
        else:
            pass
    else:
        table.type = table_type

    if tsv_metadata_formatter is not None:
        obs_md_fmt_f = observation_metadata_formatters[tsv_metadata_formatter]

    if sample_metadata is not None:
        table.add_metadata(sample_metadata)

    # if the user does not specify a name for the output metadata column,
    # set it to the same as the header key
    output_metadata_id = output_metadata_id or header_key

    if process_obs_metadata is not None and not to_tsv:
        if table.metadata(axis='observation') is None:
            raise ValueError("Observation metadata processing requested "
                             "but it doesn't appear that there is any "
                             "metadata to operate on!")

        # and if this came in as TSV, then we expect only a single type of
        # metadata
        md_key = list(table.metadata(axis='observation')[0].keys())[0]

        process_f = observation_metadata_types[process_obs_metadata]
        it = zip(table.ids(axis='observation'),
                 table.metadata(axis='observation'))
        new_md = {id_: {md_key: process_f(md[md_key])} for id_, md in it}

        if observation_metadata:
            for k, v in observation_metadata.items():
                new_md[k].update(v)
        table.add_metadata(new_md, 'observation')

    if to_tsv:
        result = table.to_tsv(header_key=header_key,
                              header_value=output_metadata_id,
                              metadata_formatter=obs_md_fmt_f)
        with open(output_filepath, 'w') as f:
            f.write(result)
        return
    elif to_json:
        fmt = 'json'
        result = table
    elif to_hdf5:
        fmt = 'hdf5'
        result = table
        if collapsed_observations:
            metadata = [{'collapsed_ids': sorted(md.keys())}
                        for md in result.metadata(axis='observation')]
            result._observation_metadata = metadata
        if collapsed_samples:
            metadata = [{'collapsed_ids': sorted(md.keys())}
                        for md in result.metadata()]
            result._sample_metadata = metadata
        if collapsed_observations or collapsed_samples:
            # We have changed the metadata, it is safer to make sure that
            # it is correct
            result._cast_metadata()
    write_biom_table(result, fmt, output_filepath)

    return
