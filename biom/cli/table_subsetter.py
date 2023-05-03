# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


import click
import h5py

from biom.cli import cli
from biom.parse import (get_axis_indices, direct_slice_data, direct_parse_key,
                        generatedby)
from biom.table import Table
from biom.util import biom_open


@cli.command(name='subset-table')
@click.option('-i', '--input-hdf5-fp', default=None,
              type=click.Path(exists=True, dir_okay=False),
              help='the input hdf5 BIOM table filepath to subset')
@click.option('-j', '--input-json-fp', default=None,
              type=click.Path(exists=True, dir_okay=False),
              help='the input json BIOM table filepath to subset')
@click.option('-a', '--axis', required=True,
              type=click.Choice(['sample', 'observation']),
              help='the axis to subset over, either sample or observation')
@click.option('-s', '--ids', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='a file containing a single column of IDs to retain '
                   '(either sample IDs or observation IDs, depending on the '
                   'axis)')
@click.option('-o', '--output-fp', required=True,
              type=click.Path(writable=True, dir_okay=False),
              help='the output BIOM table filepath')
def subset_table(input_hdf5_fp, input_json_fp, axis, ids, output_fp):
    """Subset a BIOM table.

    Subset a BIOM table, over either observations or samples, without fully
    parsing it. This command is intended to assist in working with very large
    tables when tight on memory, or as a lightweight way to subset a full
    table. Currently, it is possible to produce tables with rows or columns
    (observations or samples) that are fully zeroed.

    Example usage:

    Choose a subset of the observations in table.biom (JSON) and write them to
    subset.biom:

    $ biom subset-table -j table.biom -a observations -s observation_ids.txt \
           -o subset.biom

    Choose a subset of the observations in table.biom (HDF5) and write them to
    subset.biom:

    $ biom subset-table -i table.biom -a observations -s observation_ids.txt \
           -o subset.biom

    """
    if input_json_fp is not None:
        with open(input_json_fp) as f:
            input_json_fp = f.read()

    with open(ids) as f:
        ids = []
        for line in f:
            if not line.startswith('#'):
                ids.append(line.strip().split('\t')[0])

    table, format_ = _subset_table(input_hdf5_fp, input_json_fp, axis, ids)

    if format_ == 'json':
        with open(output_fp, 'w') as f:
            for line in table:
                f.write(line)
                f.write('\n')
    else:
        with h5py.File(output_fp, 'w') as f:
            table.to_hdf5(f, generatedby())


def _subset_table(hdf5_biom, json_table_str, axis, ids):
    if axis not in ['sample', 'observation']:
        raise ValueError("Invalid axis '%s'. Must be either 'sample' or "
                         "'observation'." % axis)

    if hdf5_biom is None and json_table_str is None:
        raise ValueError("Must specify an input table")
    elif hdf5_biom is not None and json_table_str is not None:
        raise ValueError("Can only specify one input table")

    if json_table_str is not None:
        idxs, new_axis_md = get_axis_indices(json_table_str, ids, axis)
        new_data = direct_slice_data(json_table_str, idxs, axis)

        # multiple walks over the string. bad form, but easy right now
        # ...should add a yield_and_ignore parser or something.
        def subset_generator():
            yield "{"
            yield direct_parse_key(json_table_str, "id")
            yield ","
            yield direct_parse_key(json_table_str, "format")
            yield ","
            yield direct_parse_key(json_table_str, "format_url")
            yield ","
            yield direct_parse_key(json_table_str, "type")
            yield ","
            yield direct_parse_key(json_table_str, "generated_by")
            yield ","
            yield direct_parse_key(json_table_str, "date")
            yield ","
            yield direct_parse_key(json_table_str, "matrix_type")
            yield ","
            yield direct_parse_key(json_table_str, "matrix_element_type")
            yield ","
            yield new_data
            yield ","
            yield new_axis_md
            yield ","

            if axis == "observation":
                yield direct_parse_key(json_table_str, "columns")
            else:
                yield direct_parse_key(json_table_str, "rows")
            yield "}"

        format_ = 'json'
        table = subset_generator()
    else:
        with biom_open(hdf5_biom) as f:
            table = Table.from_hdf5(f, ids=ids, axis=axis)
        format_ = 'hdf5'

    return table, format_
