#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Antonio Gonzalez",
               "Yoshiki Vazquez Baeza", "Jose Navas", "Adam Robbins-Pianka",
               "Rob Knight", "Joshua Shorenstein", "Emily TerAvest"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__version__ = "2.0.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


from .table import Table
from .parse import parse_biom_table as parse_table

example_table = Table([[0, 1, 2], [3, 4, 5]], ['O1', 'O2'],
                      ['S1', 'S2', 'S3'],
                      [{'taxonomy': ['Bacteria', 'Firmicutes']},
                       {'taxonomy': ['Bacteria', 'Bacteroidetes']}],
                      [{'environment': 'A'},
                       {'environment': 'B'},
                       {'environment': 'A'}], input_is_dense=True)


def load_table(f):
    from biom.util import biom_open
    with biom_open(f) as fp:
        table = parse_table(fp)
    return table


__all__ = ['Table', 'example_table', 'parse_table', 'load_table']

