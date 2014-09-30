#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

__author__ = "Michael Shaffer"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Michael Shaffer"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from biom.commands.table_normalizer import TableNormalizer
from biom.parse import parse_biom_table
from unittest import TestCase, main
from biom.util import HAVE_H5PY
from pyqi.core.exception import CommandError


class TableNormalizerTests(TestCase):

    def setUp(self):
        """initialize objects for use in tests"""
        self.cmd = TableNormalizer()
        self.biom_path = 'test_data/test.biom'

    def test_correct_table_type(self):
        table = self.cmd(biom_table=self.biom_path, relative_abund=True,
                         axis="sample")['table']
        if HAVE_H5PY:
            self.assertEqual(table[1], "hdf5")
        else:
            self.assertEqual(table[1], "json")

    def test_bad_inputs(self):
        # relative_abund and pa
        with self.assertRaises(CommandError):
            self.cmd(biom_table=self.biom_path, relative_abund=True,
                     presence_absence=True, axis="sample")
        # no normalization type
        with self.assertRaises(CommandError):
            self.cmd(biom_table=self.biom_path, relative_abund=False,
                     presence_absence=False, axis="sample")
        # bad axis
        with self.assertRaises(CommandError):
            self.cmd(biom_table=self.biom_path, relative_abund=True,
                     axis="nonsense")

if __name__ == "__main__":
    main()
