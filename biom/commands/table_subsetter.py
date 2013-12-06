#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from pyqi.core.command import (Command, CommandIn, CommandOut, 
        ParameterCollection)
from pyqi.core.exception import CommandError
from biom.parse import get_axis_indices, direct_slice_data, direct_parse_key
from types import GeneratorType

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__author__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class TableSubsetter(Command):
    Axes = ['samples', 'observations']

    BriefDescription = "Subset a BIOM table"
    LongDescription = ("Subset a BIOM table, over either observations or "
                       "samples, without fully parsing it. This command is "
                       "intended to assist in working with very large tables "
                       "when tight on memory, or as a lightweight way to "
                       "subset a full table. Currently, it is possible to "
                       "produce tables with rows or columns (observations or "
                       "samples) that are fully zeroed.")

    CommandIns = ParameterCollection([
        CommandIn(Name='table_str', DataType=str,
                  Description='the input BIOM table as an unparsed string',
                  Required=True),
        CommandIn(Name='axis', DataType=str,
                  Description='the axis to subset over, either ' +
                  ' or '.join(Axes), Required=True),
        CommandIn(Name='ids', DataType=list,
                  Description='the IDs to retain (either sample IDs or '
                  'observation IDs, depending on the axis)', Required=True)
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name='subset_generator',
                   DataType=GeneratorType,
                   Description='The subset generator')
    ])

    def run(self, **kwargs):
        table_str = kwargs['table_str']
        axis = kwargs['axis']
        ids = kwargs['ids']

        if axis not in self.Axes:
            raise CommandError("Invalid axis '%s'. Must be either %s." % (axis,
                ' or '.join(map(lambda e: "'%s'" % e, self.Axes))))

        idxs, new_axis_md = get_axis_indices(table_str, ids, axis)
        new_data = direct_slice_data(table_str, idxs, axis)

        # multiple walks over the string. bad form, but easy right now
        # ...should add a yield_and_ignore parser or something.
        def subset_generator():
            yield "{"
            yield direct_parse_key(table_str, "id")
            yield ","
            yield direct_parse_key(table_str, "format")
            yield ","
            yield direct_parse_key(table_str, "format_url")
            yield ","
            yield direct_parse_key(table_str, "type")
            yield ","
            yield direct_parse_key(table_str, "generated_by")
            yield ","
            yield direct_parse_key(table_str, "date")
            yield ","
            yield direct_parse_key(table_str, "matrix_type")
            yield ","
            yield direct_parse_key(table_str, "matrix_element_type")
            yield ","
            yield new_data
            yield ","
            yield new_axis_md
            yield ","

            if axis == "observations":
                yield direct_parse_key(table_str, "columns")
            else:
                yield direct_parse_key(table_str, "rows")
            yield "}"

        return {'subset_generator': subset_generator()}

CommandConstructor = TableSubsetter
