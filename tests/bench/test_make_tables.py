#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from make_tables import n_nonzero_items, get_next_row_index, get_next_col_index


class MakeTablesTests(TestCase):

    def setUp(self):
        pass

    def test_n_nonzero_items(self):
        """Determine the number of non zero items to populate"""
        n = 10
        m = 100
        p = 0.01
        exp = 10
        obs = n_nonzero_items(n, m, p)
        self.assertEqual(obs, exp)

    def test_get_next_row_index(self):
        """Determine the next row index"""
        n = 10

        cur_row = 1
        exp_row = 2
        obs_row = get_next_row_index(n, cur_row)
        self.assertEqual(obs_row, exp_row)

        cur_row = 0
        exp_row = 1
        obs_row = get_next_row_index(n, cur_row)
        self.assertEqual(obs_row, exp_row)

        cur_row = 8
        exp_row = 9
        obs_row = get_next_row_index(n, cur_row)
        self.assertEqual(obs_row, exp_row)

        cur_row = 9  # n = 10, idx of 9 is max row
        exp_row = 7
        obs_row = get_next_row_index(n, cur_row)
        self.assertEqual(obs_row, exp_row)

        cur_row = 1
        exp_row = 0
        obs_row = get_next_row_index(n, cur_row)
        self.assertEqual(obs_row, exp_row)

    def test_get_next_col_index(self):
        """determine the next col index"""
        m = 100

        cur_col = 98
        exp_col = 99
        obs_col = get_next_col_index(m, cur_col)
        self.assertEqual(obs_col, exp_col)

        cur_col = 99
        exp_col = 98
        obs_col = get_next_col_index(m, cur_col)
        self.assertEqual(obs_col, exp_col)

        cur_col = 1
        exp_col = 0
        obs_col = get_next_col_index(m, cur_col)
        self.assertEqual(obs_col, exp_col)

        cur_col = 0
        exp_col = 1
        obs_col = get_next_col_index(m, cur_col)
        self.assertEqual(obs_col, exp_col)


if __name__ == '__main__':
    main()
