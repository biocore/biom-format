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
from biom.table import Table
from biom.util import HAVE_H5PY
from biom import load_table

__author__ = "Michael Shaffer"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Michael Shaffer"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__author__ = "Michael Shaffer"
__email__ = "michael.shaffer@ucdenver.edu"


class TableNormalizer(Command):
    Axes = ['sample', 'observation']

    BriefDescription = "Normalize a BIOM table"
    LongDescription = ("Normalize the values of a BIOM table through various "
                       "methods. Relative abundance will take the relative "
                       "abundance of each observation in terms of samples or "
                       "observations.  Presence absensece will convert "
                       "observations to 1's and 0's based on presence of the "
                       "observation"
                       )

    CommandIns = ParameterCollection([
        CommandIn(Name='biom_table', DataType=str,
                  Description='the input BIOM table'),
        CommandIn(Name='axis', DataType=str,
                  Description='the axis to subset over, either ' +
                  ' or '.join(Axes),
                  Required=False),
        CommandIn(Name='relative_abund', DataType=bool,
                  Description='normalize the table by relative abundance',
                  Required=False),
        CommandIn(Name='presence_absence', DataType=bool,
                  Description='convert table to presence/absence values',
                  Required=False)
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name='table', DataType=tuple,
                   Description='The resulting table and format')
    ])

    def run(self, **kwargs):
        biom_table = kwargs['biom_table']
        axis = kwargs['axis']
        relative_abund = kwargs['relative_abund']
        p_a = kwargs['presence_absence']

        if axis not in self.Axes:
            raise CommandError("Invalid axis '%s'. Must be either %s." % (
                axis,
                ' or '.join(map(lambda e: "'%s'" % e, self.Axes))))

        if biom_table is None:
            raise CommandError("Must specify an input table")

        if relative_abund is False and p_a is False:
            raise CommandError("Must specifiy a normalization type")
        elif relative_abund is True and p_a is True:
            raise CommandError("Must specify only one normalization type")

        table = load_table(biom_table)

        if relative_abund is True:
            table.norm(axis=axis)
        else:
            table.pa()

        if HAVE_H5PY:
            return {'table': (table, 'hdf5')}
        else:
            return {'table': (table, 'json')}

CommandConstructor = TableNormalizer
