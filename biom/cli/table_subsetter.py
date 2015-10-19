# -----------------------------------------------------------------------------
# Copyright (c) 2011-2015, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division

from biom.parse import get_axis_indices, direct_slice_data, direct_parse_key
from biom.table import Table
from biom.util import biom_open


def subset_table(hdf5_biom, json_table_str, axis, ids):
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
