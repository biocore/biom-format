#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main

import os

import biom
from biom.cli.table_normalizer import _normalize_table
from biom.exception import UnknownAxisError


class TableNormalizerTests(TestCase):

    def setUp(self):
        """initialize objects for use in tests"""
        self.cmd = _normalize_table

        cwd = os.getcwd()
        if os.path.sep in __file__:
            os.chdir(os.path.dirname(__file__))
        self.table = biom.load_table(os.path.join('test_data', 'test.json'))
        os.chdir(cwd)

    def test_bad_inputs(self):
        # relative_abund and pa
        with self.assertRaises(ValueError):
            self.cmd(self.table, relative_abund=True,
                     presence_absence=True, axis="sample")
        # no normalization type
        with self.assertRaises(ValueError):
            self.cmd(self.table, relative_abund=False,
                     presence_absence=False, axis="sample")
        # bad axis
        with self.assertRaises(UnknownAxisError):
            self.cmd(self.table, relative_abund=True,
                     axis="nonsense")


if __name__ == "__main__":
    main()
