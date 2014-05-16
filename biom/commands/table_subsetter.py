#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division
from pyqi.core.command import (Command, CommandIn, CommandOut,
                               ParameterCollection)
from pyqi.core.exception import CommandError
from biom.parse import get_axis_indices, direct_slice_data, direct_parse_key
from biom.table import Table
from biom.util import biom_open

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__author__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class TableSubsetter(Command):
    Axes = ['sample', 'observation']

    BriefDescription = "Subset a BIOM table"
    LongDescription = ("Subset a BIOM table, over either observations or "
                       "samples, without fully parsing it. This command is "
                       "intended to assist in working with very large tables "
                       "when tight on memory, or as a lightweight way to "
                       "subset a full table. Currently, it is possible to "
                       "produce tables with rows or columns (observations or "
                       "samples) that are fully zeroed.")

    CommandIns = ParameterCollection([
        CommandIn(Name='json_table_str', DataType=str,
                  Description='the input BIOM table as an unparsed json '
                              'string',
                  Required=False),
        CommandIn(Name='hdf5_table', DataType=str,
                  Description='the fp to the input BIOM table',
                  Required=False),
        CommandIn(Name='axis', DataType=str,
                  Description='the axis to subset over, either ' +
                  ' or '.join(Axes), Required=True),
        CommandIn(Name='ids', DataType=list,
                  Description='the IDs to retain (either sample IDs or '
                  'observation IDs, depending on the axis)', Required=True)
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name='subsetted_table', DataType=tuple,
                   Description='The subset generator')
    ])

    def run(self, **kwargs):
        json_table_str = kwargs['json_table_str']
        hdf5_biom = kwargs['hdf5_table']
        axis = kwargs['axis']
        ids = kwargs['ids']

        if axis not in self.Axes:
            raise CommandError("Invalid axis '%s'. Must be either %s." % (
                axis,
                ' or '.join(map(lambda e: "'%s'" % e, self.Axes))))

        if hdf5_biom is None and json_table_str is None:
            raise CommandError("Must specify an input table")
        elif hdf5_biom is not None and json_table_str is not None:
            raise CommandError("Can only specify one input table")

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

        return {'subsetted_table': (table, format_)}

CommandConstructor = TableSubsetter
