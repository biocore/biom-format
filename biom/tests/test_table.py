# ----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from json import loads
from tempfile import NamedTemporaryFile
from unittest import TestCase, main
from io import StringIO
from datetime import datetime
import warnings

import numpy.testing as npt
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, csc_matrix
import scipy.sparse
import pandas.testing as pdt
import pandas as pd
import pytest

from biom import example_table, load_table, concat
from biom.exception import (UnknownAxisError, UnknownIDError, TableException,
                            DisjointIDError)
from biom.util import unzip, H5PY_VLEN_STR
from biom.table import (Table, prefer_self, index_list, list_nparray_to_sparse,
                        list_dict_to_sparse, dict_to_sparse,
                        coo_arrays_to_sparse, list_list_to_sparse,
                        nparray_to_sparse, list_sparse_to_sparse,
                        _identify_bad_value, general_parser)
from biom.parse import parse_biom_table
from biom.err import errstate
import h5py

np.random.seed(1234)


try:
    import anndata
    anndata.__version__
    HAVE_ANNDATA = True
except ImportError:
    HAVE_ANNDATA = False

try:
    import skbio
    HAVE_SKBIO = True
except ImportError:
    HAVE_SKBIO = False

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2017, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Justin Kuczynski",
               "Greg Caporaso", "Jose Clemente", "Adam Robbins-Pianka",
               "Joshua Shorenstein", "Jose Antonio Navas Molina",
               "Jorge Canardo Alastuey", "Steven Brown"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class SupportTests(TestCase):

    def test_head(self):
        # example table is 2 x 3, so no change in contained data
        exp = example_table
        obs = example_table.head()
        self.assertIsNot(obs, exp)
        self.assertEqual(obs, exp)

    def test_head_bounded(self):
        obs = example_table.head(1)
        exp = Table(np.array([[0., 1., 2.]]), ['O1'], ['S1', 'S2', 'S3'],
                    [{'taxonomy': ['Bacteria', 'Firmicutes']}],
                    [{'environment': 'A'}, {'environment': 'B'},
                     {'environment': 'A'}])

        self.assertEqual(obs, exp)

        obs = example_table.head(m=2)
        exp = Table(np.array([[0., 1.], [3., 4.]]), ['O1', 'O2'], ['S1', 'S2'],
                    [{'taxonomy': ['Bacteria', 'Firmicutes']},
                     {'taxonomy': ['Bacteria', 'Bacteroidetes']}],
                    [{'environment': 'A'}, {'environment': 'B'}])
        self.assertEqual(obs, exp)

    def test_head_overstep(self):
        # silently works
        exp = example_table
        obs = example_table.head(10000)
        self.assertIsNot(obs, exp)
        self.assertEqual(obs, exp)

    def test_head_zero_or_neg(self):
        with self.assertRaises(IndexError):
            example_table.head(0)

        with self.assertRaises(IndexError):
            example_table.head(-1)

        with self.assertRaises(IndexError):
            example_table.head(m=0)

        with self.assertRaises(IndexError):
            example_table.head(m=-1)

        with self.assertRaises(IndexError):
            example_table.head(0, 5)

        with self.assertRaises(IndexError):
            example_table.head(5, 0)

    def test_remove_empty_sample(self):
        wrn = "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient."  # noqa
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=wrn)
            t = example_table.copy()
            t._data[:, 0] = 0
            t.remove_empty()
            exp = example_table.filter({'S2', 'S3'}, inplace=False)
            self.assertEqual(t, exp)

    def test_remove_empty_obs(self):
        wrn = "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient."  # noqa
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=wrn)
            t = example_table.copy()
            t._data[0, :] = 0
            t.remove_empty()
            exp = example_table.filter({'O2', }, axis='observation',
                                       inplace=False)
            self.assertEqual(t, exp)

    def test_remove_empty_both(self):
        wrn = "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient."  # noqa
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=wrn)
            t = example_table.copy()
            t._data[:, 0] = 0
            t._data[0, :] = 0
            obs_base = t.copy()

            obs = obs_base.remove_empty(inplace=False)
            exp = example_table.filter({'S2', 'S3'}, inplace=False)
            exp = exp.filter({'O2', }, axis='observation', inplace=False)
            self.assertEqual(obs, exp)

    def test_remove_empty_identity(self):
        obs = example_table.copy()
        obs.remove_empty()
        exp = example_table.copy()
        self.assertEqual(obs, exp)

    def test_concat_table_type(self):
        table1 = example_table.copy()
        table1.type = 'foo'
        table2 = example_table.copy()
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})

        exp = Table(np.array([[0, 1, 2, 0, 1, 2],
                              [3, 4, 5, 3, 4, 5]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
                    example_table.metadata(axis='observation'),
                    list(example_table.metadata()) * 2,
                    type='foo')
        obs = table1.concat([table2, ], axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_single_table_nonlist(self):
        table2 = example_table.copy()
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})

        exp = Table(np.array([[0, 1, 2, 0, 1, 2],
                              [3, 4, 5, 3, 4, 5]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
                    example_table.metadata(axis='observation'),
                    list(example_table.metadata()) * 2)
        obs = example_table.concat(table2, axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_empty(self):
        exp = example_table.copy()
        obs = example_table.concat([])
        self.assertEqual(obs, exp)

    def test_concat_wrapper(self):
        table2 = example_table.copy()
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})

        exp = Table(np.array([[0, 1, 2, 0, 1, 2],
                              [3, 4, 5, 3, 4, 5]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
                    example_table.metadata(axis='observation'),
                    list(example_table.metadata()) * 2)
        obs = concat([example_table, table2], axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_samples(self):
        table2 = example_table.copy()
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})

        exp = Table(np.array([[0, 1, 2, 0, 1, 2],
                              [3, 4, 5, 3, 4, 5]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
                    example_table.metadata(axis='observation'),
                    list(example_table.metadata()) * 2)
        obs = example_table.concat([table2, ], axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_observations(self):
        table2 = example_table.copy()
        table2.update_ids({'O1': 'O3', 'O2': 'O4'}, axis='observation')

        exp = Table(np.array([[0, 1, 2],
                              [3, 4, 5],
                              [0, 1, 2],
                              [3, 4, 5]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3'],
                    list(example_table.metadata(axis='observation')) * 2,
                    example_table.metadata())
        obs = example_table.concat([table2, ], axis='observation')
        self.assertEqual(obs, exp)

    def test_concat_multiple(self):
        table2 = example_table.copy()
        table2.update_ids({'O1': 'O3', 'O2': 'O4'}, axis='observation')
        table3 = example_table.copy()
        table3.update_ids({'O1': 'O5', 'O2': 'O6'}, axis='observation')

        exp = Table(np.array([[0, 1, 2],
                              [3, 4, 5],
                              [0, 1, 2],
                              [3, 4, 5],
                              [0, 1, 2],
                              [3, 4, 5]]),
                    ['O1', 'O2', 'O3', 'O4', 'O5', 'O6'],
                    ['S1', 'S2', 'S3'],
                    list(example_table.metadata(axis='observation')) * 3,
                    example_table.metadata())
        obs = example_table.concat([table2, table3], axis='observation')
        self.assertEqual(obs, exp)

    def test_concat_different_order(self):
        table2 = example_table.sort_order(['S3', 'S2', 'S1'])
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})
        table2 = table2.sort_order(['O2', 'O1'], axis='observation')
        table3 = example_table.sort_order(['S2', 'S1', 'S3'])
        table3.update_ids({'S1': 'S7', 'S2': 'S8', 'S3': 'S9'})

        exp = Table(np.array([[0, 1, 2, 2, 1, 0, 1, 0, 2],
                              [3, 4, 5, 5, 4, 3, 4, 3, 5]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3', 'S6', 'S5', 'S4', 'S8', 'S7', 'S9'],
                    example_table.metadata(axis='observation'),
                    list(example_table.metadata()) +
                    list(table2.metadata()) +
                    list(table3.metadata()))
        obs = example_table.concat([table2, table3], axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_pad_on_subset(self):
        table2 = example_table.copy()
        table2.update_ids({'O2': 'O3'}, axis='observation', strict=False)
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})
        exp_obs_md = list(example_table.metadata(axis='observation'))
        exp_obs_md.append(example_table.metadata('O2', axis='observation'))

        exp = Table(np.array([[0, 1, 2, 0, 1, 2],
                              [3, 4, 5, 0, 0, 0],
                              [0, 0, 0, 3, 4, 5]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
                    exp_obs_md,
                    list(example_table.metadata()) * 2)

        obs = example_table.concat([table2, ], axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_no_metadata_bug(self):
        table1 = example_table.copy()
        table1._sample_metadata = None
        table1._observation_metadata = None
        table2 = example_table.copy()
        table2._sample_metadata = None
        table2._observation_metadata = None
        table2.update_ids({'O2': 'O3'}, axis='observation', strict=False)
        table2.update_ids({'S1': 'S4', 'S2': 'S5', 'S3': 'S6'})

        exp = Table(np.array([[0, 1, 2, 0, 1, 2],
                              [3, 4, 5, 0, 0, 0],
                              [0, 0, 0, 3, 4, 5]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])

        obs = table1.concat([table2, ], axis='sample')
        self.assertEqual(obs, exp)

    def test_concat_raise_overlap(self):
        with self.assertRaises(DisjointIDError):
            example_table.concat([example_table])

        with self.assertRaises(DisjointIDError):
            example_table.concat([example_table], axis='observation')

    def test_align_to_no_overlap(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['w', 'x'], ['y', 'z'])

        with self.assertRaises(DisjointIDError):
            a.align_to(b)

    def test_align_to_overlap_observation(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['b', 'a'], ['y', 'z'])
        exp = Table(np.array([[2, 3], [0, 1]]), ['a', 'b'], ['y', 'z'])
        obs = b.align_to(a, axis='observation')
        self.assertEqual(obs, exp)

    def test_align_to_overlap_observation_no_overlap(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['x', 'y'], ['y', 'z'])
        with self.assertRaises(DisjointIDError):
            b.align_to(a, axis='observation')

    def test_align_to_overlap_sample_no_overlap(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['b', 'a'], ['y', 'z'])
        with self.assertRaises(DisjointIDError):
            b.align_to(a, axis='sample')

    def test_align_to_overlap_both_no_overlap(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['b', 'a'], ['y', 'z'])
        with self.assertRaises(DisjointIDError):
            # should fail if one axis doesn't overlap
            b.align_to(a, axis='both')

    def test_align_to_overlap_autodetect_observation(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['b', 'a'], ['y', 'z'])
        exp = Table(np.array([[2, 3], [0, 1]]), ['a', 'b'], ['y', 'z'])
        obs = b.align_to(a)
        self.assertEqual(obs, exp)

    def test_align_to_overlap_sample(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['d', 'c'])
        exp = Table(np.array([[1, 0], [3, 2]]), ['a', 'b'], ['c', 'd'])
        obs = b.align_to(a, axis='sample')
        self.assertEqual(obs, exp)

    def test_align_to_overlap_autodetect_sample(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['d', 'c'])
        exp = Table(np.array([[1, 0], [3, 2]]), ['a', 'b'], ['c', 'd'])
        obs = b.align_to(a)
        self.assertEqual(obs, exp)

    def test_align_to_overlap_both(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['b', 'a'], ['d', 'c'])
        exp = Table(np.array([[3, 2], [1, 0]]), ['a', 'b'], ['c', 'd'])
        obs = b.align_to(a, axis='both')
        self.assertEqual(obs, exp)

    def test_align_to_overlap_autodetect_both(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[0, 1], [2, 3]]), ['b', 'a'], ['d', 'c'])
        exp = Table(np.array([[3, 2], [1, 0]]), ['a', 'b'], ['c', 'd'])
        obs = b.align_to(a)
        self.assertEqual(obs, exp)

    def test_align_to_overlap_autodetect_vary_values(self):
        a = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'])
        b = Table(np.array([[10, 11], [12, 13]]), ['b', 'a'], ['d', 'c'])
        exp = Table(np.array([[13, 12], [11, 10]]), ['a', 'b'], ['c', 'd'])
        obs = b.align_to(a)
        self.assertEqual(obs, exp)

    def test_table_sparse_nparray(self):
        """beat the table sparsely to death"""
        # nparray test
        samp_ids = ['1', '2', '3', '4']
        obs_ids = ['a', 'b', 'c']
        nparray = np.array([[1, 2, 3, 4], [-1, 6, 7, 8], [9, 10, 11, 12]])
        data = nparray_to_sparse(
            np.array([[1, 2, 3, 4], [-1, 6, 7, 8], [9, 10, 11, 12]]))
        exp = Table(data, obs_ids, samp_ids)
        obs = Table(nparray, obs_ids, samp_ids)
        self.assertEqual(obs, exp)

    def test_table_sparse_list_nparray(self):
        """beat the table sparsely to death"""
        # list of nparray test
        samp_ids = ['1', '2', '3', '4']
        obs_ids = ['a', 'b', 'c']
        list_np = [np.array([1, 2, 3, 4]), np.array([5, 6, 7, 8]),
                   np.array([9, 10, 11, 12])]
        data = list_nparray_to_sparse(list_np)
        exp = Table(data, obs_ids, samp_ids)
        obs = Table(list_np, obs_ids, samp_ids)
        self.assertEqual(obs, exp)

    def test_table_sparse_dict(self):
        """beat the table sparsely to death"""
        # dict test
        samp_ids = range(24)
        obs_ids = range(101)
        dict_input = {(0, 0): 1, (0, 10): 5, (100, 23): -3}
        d_input = np.zeros((101, 24), dtype=float)
        d_input[0, 0] = 1
        d_input[0, 10] = 5
        d_input[100, 23] = -3
        data = nparray_to_sparse(d_input)
        exp = Table(data, obs_ids, samp_ids)
        obs = Table(dict_input, obs_ids, samp_ids)
        self.assertEqual(obs, exp)

    def test_table_sparse_list_dict(self):
        """beat the table sparsely to death"""
        # list of dict test
        samp_ids = range(11)
        obs_ids = range(3)
        ld_input = np.zeros((3, 11), dtype=float)
        ld_input[0, 5] = 10
        ld_input[0, 10] = 2
        ld_input[1, 1] = 15
        ld_input[2, 3] = 7
        data = nparray_to_sparse(ld_input)
        exp = Table(data, obs_ids, samp_ids)
        list_dict = [{(0, 5): 10, (10, 10): 2}, {(0, 1): 15}, {(0, 3): 7}]
        obs = Table(list_dict, obs_ids, samp_ids)
        self.assertEqual(obs, exp)

    def test_table_sparse_list_list(self):
        """beat the table sparsely to death"""
        # list list test
        samp_ids = range(3)
        obs_ids = range(2)
        exp_data = lil_matrix((2, 3))
        exp_data[0, 1] = 5
        exp_data[1, 2] = 10
        exp = Table(exp_data, obs_ids, samp_ids)
        input_ = [[0, 1, 5], [1, 2, 10]]
        obs = Table(input_, obs_ids, samp_ids)
        self.assertEqual(obs, exp)

    def test_table_exception(self):
        """Make sure a TableException can be raised"""
        def f():
            raise TableException
        self.assertRaises(TableException, f)

    def test_prefer_self(self):
        """prefer x"""
        exp = 1
        obs = prefer_self(1, 2)
        self.assertEqual(obs, exp)

        exp = 2
        obs = prefer_self(None, 2)
        self.assertEqual(obs, exp)

        exp = None
        obs = prefer_self(None, None)
        self.assertEqual(obs, exp)

    def test_index_list(self):
        """returns a dict for list lookups"""
        exp = {'a': 2, 'b': 0, 'c': 1}
        obs = index_list(['b', 'c', 'a'])
        self.assertEqual(obs, exp)


class TableTests(TestCase):

    def setUp(self):
        self.simple_derived = Table(
            np.array([[5, 6], [7, 8]]), [3, 4], [1, 2])
        self.vals = {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8}
        self.st1 = Table(self.vals, ['1', '2'], ['a', 'b'])
        self.st2 = Table(self.vals, ['1', '2'], ['a', 'b'])
        self.vals3 = {(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}
        self.vals4 = {(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}
        self.st3 = Table(self.vals3, ['2', '3'], ['b', 'c'])
        self.st4 = Table(self.vals4, ['3', '4'],  ['c', 'd'])
        self.st_rich = Table(self.vals,
                             ['1', '2'], ['a', 'b'],
                             [{'taxonomy': ['k__a', 'p__b']},
                              {'taxonomy': ['k__a', 'p__c']}],
                             [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                             )
        self.st_group_rich = Table(
            self.vals, ['1', '2'], ['a', 'b'],
            [{'taxonomy': ['k__a', 'p__b']}, {'taxonomy': ['k__a', 'p__c']}],
            [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
            observation_group_metadata={'tree': ('newick', '(a:0.3,b:0.4);')},
            sample_group_metadata={'category': ('newick', '(1:0.3,2:0.4);')}
            )

        self.empty_st = Table([], [], [])

        self.vals5 = {(0, 1): 2, (1, 1): 4}
        self.st5 = Table(self.vals5, ['5', '6'], ['a', 'b'])

        self.vals6 = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
        self.st6 = Table(self.vals6, ['5', '6'], ['a', 'b'])

        self.vals7 = {(0, 0): 5, (0, 1): 7, (1, 0): 8, (1, 1): 0}
        self.st7 = Table(self.vals7, ['5', '6'], ['a', 'b'])

        self.single_sample_st = Table(
            np.array([[2.0], [0.0], [1.0]]), ['O1', 'O2', 'O3'],
            ['S1'])
        self.single_obs_st = Table(np.array([[2.0, 0.0, 1.0]]),
                                   ['01'], ['S1', 'S2', 'S3'])

        self.to_remove = []

        # 1 0 2
        # 3 0 4
        self.mat1 = Table(np.array([[1, 0, 2], [3, 0, 4]]),
                          ['o1', 'o2'], ['s1', 's2', 's3'])

        # Empty/null cases (i.e., 0x0, 0xn, nx0).
        def ids(X):
            return ['x%d' % e for e in range(0, X)]
        self.null1 = Table(np.zeros((0, 0)), [], [])
        self.null2 = Table(
            np.zeros((0, 42), dtype=float), [], ids(42))
        self.null3 = Table(
            np.zeros((42, 0), dtype=float), ids(42), [])
        self.nulls = [self.null1, self.null2, self.null3]

        # 0 0
        # 0 0
        self.empty = Table(np.zeros((2, 2)), ids(2), ids(2))

        # 1 0 3
        h = np.array([[1.0, 0.0, 3.0]])
        self.row_vec = Table(h, ids(1), ids(3))

        # 1
        # 0
        # 3
        h = np.array([[1], [0], [3]])
        self.col_vec = Table(h, ids(3), ids(1))

        # 1x1
        h = np.array([[42]])
        self.single_ele = Table(h, ['b'], ['a'])

        # Explicit zeros.
        self.explicit_zeros = Table(np.array([[0, 0, 1], [1, 0, 0],
                                              [1, 0, 2]]),
                                    ['a', 'b', 'c'], ['x', 'y', 'z'])

    def tearDown(self):
        if self.to_remove:
            for f in self.to_remove:
                os.remove(f)

    def test_data_property(self):
        exp = self.simple_derived._data
        obs = self.simple_derived.matrix_data
        self.assertEqual((obs != exp).nnz, 0)

        with self.assertRaises(AttributeError):
            self.simple_derived.matrix_data = 'foo'

    def test_repr(self):
        """__repr__ method of biom.table.Table"""
        # table
        data = np.asarray([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
        t = Table(data, ['a', 'b', 'c'], ['x', 'y', 'z'])
        self.assertEqual("3 x 3 <class 'biom.table.Table'> with 1 nonzero "
                         "entries (11% dense)", repr(t))

        # empty table
        data = np.asarray([[]])
        t = Table(data, [], [])
        self.assertEqual("0 x 0 <class 'biom.table.Table'> with 0 nonzero "
                         "entries (0% dense)", repr(t))

    def test_init_with_nparray(self):
        """to_sparse in constructor should be triggered"""
        data = np.array([[1, 2], [3, 4]])
        samp_ids = ['a', 'b']
        obs_ids = ['1', '2']
        exp = Table(data, obs_ids, samp_ids)
        obs = Table(data, obs_ids, samp_ids)
        self.assertEqual(obs, exp)

    def test_min_observation(self):
        exp = np.array([5, 7])
        obs = self.simple_derived.min('observation')
        npt.assert_equal(obs, exp)

    def test_min_sample(self):
        exp = np.array([5, 6])
        obs = self.simple_derived.min('sample')
        npt.assert_equal(obs, exp)

    def test_min_whole(self):
        exp = 5
        obs = self.simple_derived.min('whole')
        npt.assert_equal(obs, exp)

    def test_max_observation(self):
        exp = np.array([6, 8])
        obs = self.simple_derived.max('observation')
        npt.assert_equal(obs, exp)

    def test_max_sample(self):
        exp = np.array([7, 8])
        obs = self.simple_derived.max('sample')
        npt.assert_equal(obs, exp)

    def test_max_whole(self):
        exp = 8
        obs = self.simple_derived.max('whole')
        npt.assert_equal(obs, exp)

    def test_general_parser(self):
        test_and_exp = [(b'foo', 'foo'),
                        ('foo', 'foo'),
                        (b'', ''),
                        ('', ''),
                        (b'10', '10'),
                        ('10', '10'),
                        (b'3.14', '3.14'),
                        ('3.14', '3.14')]
        for test, exp in test_and_exp:
            obs = general_parser(test)
            self.assertEqual(obs, exp)

    def test_from_hdf5_non_hdf5_file_or_group(self):
        with self.assertRaises(ValueError):
            Table.from_hdf5(10)

    def test_from_hdf5_empty_md(self):
        """Parse a hdf5 formatted BIOM table w/o metadata"""
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/empty.biom', 'r'))
        os.chdir(cwd)

        self.assertTrue(t._sample_metadata is None)
        self.assertTrue(t._observation_metadata is None)

    def test_from_hdf5_custom_parsers(self):
        def parser(item):
            return general_parser(item).upper()
        parse_fs = {'BODY_SITE': parser}

        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'),
                            parse_fs=parse_fs)
        os.chdir(cwd)

        for m in t.metadata():
            self.assertIn(m['BODY_SITE'], ('GUT', 'SKIN', b'GUT', b'SKIN'))

    def test_from_hdf5_issue_731(self):
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'))
        os.chdir(cwd)
        self.assertTrue(isinstance(t.table_id, str))
        self.assertTrue(isinstance(t.type, str))

    def test_from_hdf5_group_metadata_issue_926(self):
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test_grp_metadata.biom', 'r'))
        os.chdir(cwd)
        obs = t.group_metadata(axis='observation')
        self.assertEqual(obs, {'tree': '(GG_OTU_5:0.1,(GG_OTU_3:0.1,GG_OTU_4:0.5):0.1,(GG_OTU_1:0.3,GG_OTU_2:0.4):0.1);'})  # noqa

    def test_from_hdf5(self):
        """Parse a hdf5 formatted BIOM table"""
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'))
        os.chdir(cwd)

        npt.assert_equal(t.ids(), ('Sample1', 'Sample2', 'Sample3',
                                   'Sample4', 'Sample5', 'Sample6'))
        npt.assert_equal(t.ids(axis='observation'),
                         ('GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3',
                          'GG_OTU_4', 'GG_OTU_5'))
        exp_obs_md = (
            {'taxonomy': [
                'k__Bacteria',
                'p__Proteobacteria',
                'c__Gammaproteobacteria',
                'o__Enterobacteriales',
                'f__Enterobacteriaceae',
                'g__Escherichia',
                's__',
            ]},
            {'taxonomy': [
                'k__Bacteria',
                'p__Cyanobacteria',
                'c__Nostocophycideae',
                'o__Nostocales',
                'f__Nostocaceae',
                'g__Dolichospermum',
                's__',
            ]},
            {'taxonomy': [
                'k__Archaea',
                'p__Euryarchaeota',
                'c__Methanomicrobia',
                'o__Methanosarcinales',
                'f__Methanosarcinaceae',
                'g__Methanosarcina',
                's__',
            ]},
            {'taxonomy': [
                'k__Bacteria',
                'p__Firmicutes',
                'c__Clostridia',
                'o__Halanaerobiales',
                'f__Halanaerobiaceae',
                'g__Halanaerobium',
                's__Halanaerobiumsaccharolyticum',
            ]},
            {'taxonomy': [
                'k__Bacteria',
                'p__Proteobacteria',
                'c__Gammaproteobacteria',
                'o__Enterobacteriales',
                'f__Enterobacteriaceae',
                'g__Escherichia',
                's__',
            ]})
        self.assertEqual(t._observation_metadata, exp_obs_md)

        exp_samp_md = ({'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CGCTTATCGAGA',
                        'Description': 'human gut',
                        'BODY_SITE': 'gut'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CATACCAGTAGC',
                        'Description': 'human gut',
                        'BODY_SITE': 'gut'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTCTCTACCTGT',
                        'Description': 'human gut',
                        'BODY_SITE': 'gut'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTCTCGGCCTGT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTCTCTACCAAT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTAACTACCAAT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'})
        self.assertEqual(t._sample_metadata, exp_samp_md)

        exp = [np.array([0., 0., 1., 0., 0., 0.]),
               np.array([5., 1., 0., 2., 3., 1.]),
               np.array([0., 0., 1., 4., 0., 2.]),
               np.array([2., 1., 1., 0., 0., 1.]),
               np.array([0., 1., 1., 0., 0., 0.])]
        npt.assert_equal(list(t.iter_data(axis="observation")), exp)

    def test_from_hdf5_sample_subset_no_metadata(self):
        """Parse a sample subset of a hdf5 formatted BIOM table"""
        samples = [b'Sample2', b'Sample4', b'Sample6']

        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'), ids=samples,
                            subset_with_metadata=False)
        os.chdir(cwd)

        npt.assert_equal(t.ids(), ['Sample2', 'Sample4', 'Sample6'])
        npt.assert_equal(t.ids(axis='observation'),
                         ['GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3', 'GG_OTU_4',
                          'GG_OTU_5'])
        exp_obs_md = None
        self.assertEqual(t._observation_metadata, exp_obs_md)
        exp_samp_md = None
        self.assertEqual(t._sample_metadata, exp_samp_md)

        exp = [np.array([0, 0, 0]),
               np.array([1., 2., 1.]),
               np.array([0., 4., 2.]),
               np.array([1., 0., 1.]),
               np.array([1., 0., 0.])]
        npt.assert_equal(list(t.iter_data(axis='observation')), exp)

    def test_from_hdf5_sample_subset(self):
        """Parse a sample subset of a hdf5 formatted BIOM table"""
        samples = ['Sample2', 'Sample4', 'Sample6']

        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'), ids=samples)
        os.chdir(cwd)

        npt.assert_equal(t.ids(), ['Sample2', 'Sample4', 'Sample6'])
        npt.assert_equal(t.ids(axis='observation'),
                         ['GG_OTU_2', 'GG_OTU_3', 'GG_OTU_4', 'GG_OTU_5'])
        exp_obs_md = (
            {'taxonomy': [
                'k__Bacteria',
                'p__Cyanobacteria',
                'c__Nostocophycideae',
                'o__Nostocales',
                'f__Nostocaceae',
                'g__Dolichospermum',
                's__',
            ]},
            {'taxonomy': [
                'k__Archaea',
                'p__Euryarchaeota',
                'c__Methanomicrobia',
                'o__Methanosarcinales',
                'f__Methanosarcinaceae',
                'g__Methanosarcina',
                's__',
            ]},
            {'taxonomy': [
                'k__Bacteria',
                'p__Firmicutes',
                'c__Clostridia',
                'o__Halanaerobiales',
                'f__Halanaerobiaceae',
                'g__Halanaerobium',
                's__Halanaerobiumsaccharolyticum',
            ]},
            {'taxonomy': [
                'k__Bacteria',
                'p__Proteobacteria',
                'c__Gammaproteobacteria',
                'o__Enterobacteriales',
                'f__Enterobacteriaceae',
                'g__Escherichia',
                's__',
            ]})
        self.assertEqual(t._observation_metadata, exp_obs_md)

        exp_samp_md = ({'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CATACCAGTAGC',
                        'Description': 'human gut',
                        'BODY_SITE': 'gut'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTCTCGGCCTGT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTAACTACCAAT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'})
        self.assertEqual(t._sample_metadata, exp_samp_md)

        exp = [np.array([1., 2., 1.]),
               np.array([0., 4., 2.]),
               np.array([1., 0., 1.]),
               np.array([1., 0., 0.])]
        npt.assert_equal(list(t.iter_data(axis='observation')), exp)

    def test_from_hdf5_observation_subset_no_metadata(self):
        """Parse a observation subset of a hdf5 formatted BIOM table"""
        observations = [b'GG_OTU_1', b'GG_OTU_3', b'GG_OTU_5']

        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'),
                            ids=observations, axis='observation',
                            subset_with_metadata=False)
        os.chdir(cwd)

        npt.assert_equal(t.ids(), ['Sample1', 'Sample2', 'Sample3',
                                   'Sample4', 'Sample5', 'Sample6'])
        npt.assert_equal(t.ids(axis='observation'),
                         ['GG_OTU_1', 'GG_OTU_3', 'GG_OTU_5'])
        exp_obs_md = None
        self.assertEqual(t._observation_metadata, exp_obs_md)

        exp_samp_md = None
        self.assertEqual(t._sample_metadata, exp_samp_md)

        exp = [np.array([0, 0., 1., 0., 0, 0.]),
               np.array([0, 0., 1., 4., 0, 2.]),
               np.array([0, 1., 1., 0., 0, 0.])]
        npt.assert_equal(list(t.iter_data(axis='observation')), exp)

    def test_from_hdf5_observation_subset(self):
        """Parse a observation subset of a hdf5 formatted BIOM table"""
        observations = ['GG_OTU_1', 'GG_OTU_3', 'GG_OTU_5']

        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/test.biom', 'r'),
                            ids=observations, axis='observation')
        os.chdir(cwd)

        npt.assert_equal(t.ids(), ['Sample2', 'Sample3', 'Sample4',
                                   'Sample6'])
        npt.assert_equal(t.ids(axis='observation'),
                         ['GG_OTU_1', 'GG_OTU_3', 'GG_OTU_5'])
        exp_obs_md = (
            {'taxonomy': [
                'k__Bacteria',
                'p__Proteobacteria',
                'c__Gammaproteobacteria',
                'o__Enterobacteriales',
                'f__Enterobacteriaceae',
                'g__Escherichia',
                's__',
            ]},
            {'taxonomy': [
                'k__Archaea',
                'p__Euryarchaeota',
                'c__Methanomicrobia',
                'o__Methanosarcinales',
                'f__Methanosarcinaceae',
                'g__Methanosarcina',
                's__',
            ]},
            {'taxonomy': [
                'k__Bacteria',
                'p__Proteobacteria',
                'c__Gammaproteobacteria',
                'o__Enterobacteriales',
                'f__Enterobacteriaceae',
                'g__Escherichia',
                's__',
            ]})
        self.assertEqual(t._observation_metadata, exp_obs_md)

        exp_samp_md = ({'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CATACCAGTAGC',
                        'Description': 'human gut',
                        'BODY_SITE': 'gut'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTCTCTACCTGT',
                        'Description': 'human gut',
                        'BODY_SITE': 'gut'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTCTCGGCCTGT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'},
                       {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                        'BarcodeSequence': 'CTAACTACCAAT',
                        'Description': 'human skin',
                        'BODY_SITE': 'skin'})
        self.assertEqual(t._sample_metadata, exp_samp_md)

        exp = [np.array([0., 1., 0., 0.]),
               np.array([0., 1., 4., 2.]),
               np.array([1., 1., 0., 0.])]
        npt.assert_equal(list(t.iter_data(axis='observation')), exp)

    def test_from_hdf5_subset_error(self):
        """hdf5 biom table parse throws error with invalid parameters"""
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])

        # Raises an error if not all the given samples are in the OTU table
        with self.assertRaises(ValueError):
            Table.from_hdf5(h5py.File('test_data/test.biom', 'r'),
                            ids=['Sample2', 'DoesNotExist', 'Sample6'])

        # Raises an error if not all the given observation are in the OTU table
        with self.assertRaises(ValueError):
            Table.from_hdf5(h5py.File('test_data/test.biom', 'r'),
                            ids=['GG_OTU_1', 'DoesNotExist'],
                            axis='observation')
        os.chdir(cwd)

    def test_from_hdf5_empty_table(self):
        """HDF5 biom parse successfully loads an empty table"""
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        t = Table.from_hdf5(h5py.File('test_data/empty.biom', 'r'))
        os.chdir(cwd)

        npt.assert_equal(t.ids(), [])
        npt.assert_equal(t.ids(axis='observation'), [])
        self.assertEqual(t._observation_metadata, None)
        self.assertEqual(t._sample_metadata, None)
        npt.assert_equal(list(t.iter_data(axis='observation')), [])

    def test_to_from_hdf5_bug_861(self):
        t = Table(np.array([[0, 1, 2], [3, 4, 5]]), ['a', 'b'],
                  ['c', 'd', 'e'])
        t.add_metadata({'a': {'a / problem': 10}, 'b': {'a / problem': 20}},
                       axis='observation')
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            t.to_hdf5(h5, 'tests')
            h5.close()

            h5 = h5py.File(tmpfile.name, 'r')
            obs = Table.from_hdf5(h5)
            h5.close()
        self.to_remove.append(tmpfile.name)

        self.assertEqual(obs, t)

    def test_to_from_hdf5_creation_date(self):
        t = Table(np.array([[0, 1, 2], [3, 4, 5]]), ['a', 'b'],
                  ['c', 'd', 'e'])
        current = datetime.now()
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            t.to_hdf5(h5, 'tests', creation_date=current)
            h5.close()

            h5 = h5py.File(tmpfile.name, 'r')
            obs = Table.from_hdf5(h5)
            self.assertEqual(obs.create_date, current)
            h5.close()
        self.to_remove.append(tmpfile.name)

        self.assertEqual(obs, t)

    def test_to_hdf5_empty_table(self):
        """Successfully writes an empty OTU table in HDF5 format"""
        # Create an empty OTU table
        t = Table([], [], [])
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            t.to_hdf5(h5, 'tests')
            h5.close()
        self.to_remove.append(tmpfile.name)

    def test_to_hdf5_empty_table_bug_619(self):
        """Successfully writes an empty OTU table in HDF5 format"""
        t = example_table.filter({}, axis='observation', inplace=False)
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            t.to_hdf5(h5, 'tests')
            h5.close()
        self.to_remove.append(tmpfile.name)

        t = example_table.filter({}, inplace=False)
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            t.to_hdf5(h5, 'tests')
            h5.close()
        self.to_remove.append(tmpfile.name)

    def test_to_hdf5_missing_metadata_observation(self):
        # exercises a vlen_list
        t = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'],
                  [{'taxonomy': None},
                   {'taxonomy': ['foo', 'baz']}])

        with NamedTemporaryFile(delete=False) as tmpfile:
            with h5py.File(tmpfile.name, 'w') as h5:
                t.to_hdf5(h5, 'tests')
            obs = load_table(tmpfile.name)
        self.to_remove.append(tmpfile.name)
        self.assertEqual(obs.metadata(axis='observation'),
                         ({'taxonomy': None},
                          {'taxonomy': ['foo', 'baz']}))

    def test_to_hdf5_missing_metadata_sample(self):
        # exercises general formatter
        t = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'], None,
                  [{'dat': None},
                   {'dat': 'foo'}])

        with NamedTemporaryFile(delete=False) as tmpfile:
            with h5py.File(tmpfile.name, 'w') as h5:
                t.to_hdf5(h5, 'tests')
            obs = load_table(tmpfile.name)
        self.to_remove.append(tmpfile.name)
        self.assertEqual(obs.metadata(axis='sample'),
                         ({'dat': ''},
                          {'dat': 'foo'}))

    def test_to_hdf5_inconsistent_metadata_categories_observation(self):
        t = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'],
                  [{'taxonomy_A': 'foo; bar'},
                   {'taxonomy_B': 'foo; baz'}])

        with NamedTemporaryFile(delete=False) as tmpfile:
            with h5py.File(tmpfile.name, 'w') as h5:
                with self.assertRaisesRegex(ValueError,
                                            'inconsistent metadata'):
                    t.to_hdf5(h5, 'tests')
        self.to_remove.append(tmpfile.name)

    def test_to_hdf5_inconsistent_metadata_categories_sample(self):
        t = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'],
                  None,
                  [{'dat_A': 'foo; bar'},
                   {'dat_B': 'foo; baz'}])

        with NamedTemporaryFile(delete=False) as tmpfile:
            with h5py.File(tmpfile.name, 'w') as h5:
                with self.assertRaisesRegex(ValueError,
                                            'inconsistent metadata'):
                    t.to_hdf5(h5, 'tests')
        self.to_remove.append(tmpfile.name)

    def test_to_hdf5_malformed_taxonomy(self):
        t = Table(np.array([[0, 1], [2, 3]]), ['a', 'b'], ['c', 'd'],
                  [{'taxonomy': 'foo; bar'},
                   {'taxonomy': 'foo; baz'}])

        with NamedTemporaryFile(delete=False) as tmpfile:
            with h5py.File(tmpfile.name, 'w') as h5:
                t.to_hdf5(h5, 'tests')
            obs = load_table(tmpfile.name)
        self.to_remove.append(tmpfile.name)
        self.assertEqual(obs.metadata(axis='observation'),
                         ({'taxonomy': ['foo', 'bar']},
                          {'taxonomy': ['foo', 'baz']}))

    def test_to_hdf5_general_fallback_to_list(self):
        st_rich = Table(self.vals,
                        ['1', '2'], ['a', 'b'],
                        [{'foo': ['k__a', 'p__b']},
                         {'foo': ['k__a', 'p__c']}],
                        [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            st_rich.to_hdf5(h5, 'tests')
            h5.close()
        self.to_remove.append(tmpfile.name)

    def test_to_hdf5_custom_formatters(self):
        self.st_rich = Table(self.vals,
                             ['1', '2'], ['a', 'b'],
                             [{'taxonomy': ['k__a', 'p__b']},
                              {'taxonomy': ['k__a', 'p__c']}],
                             [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])

        def bc_formatter(grp, category, md, compression):
            name = 'metadata/%s' % category
            data = np.array([m[category].upper().encode('utf8') for m in md])
            grp.create_dataset(name, shape=data.shape, dtype=H5PY_VLEN_STR,
                               data=data, compression=compression)

        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            self.st_rich.to_hdf5(h5, 'tests',
                                 format_fs={'barcode': bc_formatter})
            h5.close()

            h5 = h5py.File(tmpfile.name, 'r')
            self.assertIn('observation', h5)
            self.assertIn('sample', h5)
            self.assertEqual(sorted(h5.attrs.keys()), sorted(['id', 'type',
                                                              'format-url',
                                                              'format-version',
                                                              'generated-by',
                                                              'creation-date',
                                                              'shape', 'nnz']))

            obs = Table.from_hdf5(h5)
            for m1, m2 in zip(obs.metadata(), self.st_rich.metadata()):
                self.assertNotEqual(m1['barcode'], m2['barcode'])
                self.assertEqual(m1['barcode'].lower(), m2['barcode'])
            h5.close()
        self.to_remove.append(tmpfile.name)

    def test_to_hdf5(self):
        """Write a file"""
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            self.st_rich.to_hdf5(h5, 'tests')
            h5.close()

            h5 = h5py.File(tmpfile.name, 'r')
            self.assertIn('observation', h5)
            self.assertIn('sample', h5)
            self.assertEqual(sorted(h5.attrs.keys()), sorted(['id', 'type',
                                                              'format-url',
                                                              'format-version',
                                                              'generated-by',
                                                              'creation-date',
                                                              'shape', 'nnz']))

            obs = Table.from_hdf5(h5)
            self.assertEqual(obs, self.st_rich)
            h5.close()
        self.to_remove.append(tmpfile.name)

        # Test with a collapsed table
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            dt_rich = Table(
                np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
                ['1', '2', '3'], ['a', 'b', 'c'],
                [{'taxonomy': ['k__a', 'p__b']},
                 {'taxonomy': ['k__a', 'p__c']},
                 {'taxonomy': ['k__a', 'p__c']}],
                [{'barcode': 'aatt'},
                 {'barcode': 'ttgg'},
                 {'barcode': 'aatt'}])

            def bin_f(id_, x):
                return x['barcode']

            collapsed = dt_rich.collapse(
                bin_f, norm=False, min_group_size=1,
                axis='sample').sort(axis='sample')
            collapsed.to_hdf5(h5, 'tests')
            h5.close()

            h5 = h5py.File(tmpfile.name, 'r')
            self.assertIn('observation', h5)
            self.assertIn('sample', h5)
            self.assertEqual(sorted(h5.attrs.keys()), sorted(['id', 'type',
                                                              'format-url',
                                                              'format-version',
                                                              'generated-by',
                                                              'creation-date',
                                                              'shape', 'nnz']))

            obs = Table.from_hdf5(h5)
            h5.close()

            exp = Table(
                np.array([[12, 6], [18, 9], [24, 12]]),
                ['1', '2', '3'], ['aatt', 'ttgg'],
                [{'taxonomy': ['k__a', 'p__b']},
                 {'taxonomy': ['k__a', 'p__c']},
                 {'taxonomy': ['k__a', 'p__c']}],
                [{'collapsed_ids': ['a', 'c']},
                 {'collapsed_ids': ['b']}])
            self.assertEqual(obs, exp)
        self.to_remove.append(tmpfile.name)

        # Test with table having a None on taxonomy
        with NamedTemporaryFile(delete=False) as tmpfile:
            h5 = h5py.File(tmpfile.name, 'w')
            t = Table(self.vals, ['1', '2'], ['a', 'b'],
                      [{'taxonomy': ['k__a', 'p__b']},
                       {'taxonomy': None}],
                      [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])
            t.to_hdf5(h5, 'tests')
            h5.close()

            h5 = h5py.File(tmpfile.name, 'r')
            self.assertIn('observation', h5)
            self.assertIn('sample', h5)
            self.assertEqual(sorted(h5.attrs.keys()), sorted(['id', 'type',
                                                              'format-url',
                                                              'format-version',
                                                              'generated-by',
                                                              'creation-date',
                                                              'shape', 'nnz']))

            obs = Table.from_hdf5(h5)
            h5.close()
            self.assertEqual(obs, t)
        self.to_remove.append(tmpfile.name)

    def test_from_tsv(self):
        tab1_fh = StringIO(otu_table1)
        sparse_rich = Table.from_tsv(tab1_fh, None, None,
                                     OBS_META_TYPES['naive'])
        self.assertEqual(sorted(sparse_rich.ids()),
                         sorted(['Fing', 'Key', 'NA']))
        self.assertEqual(sorted(sparse_rich.ids(axis='observation')),
                         list(map(str, [0, 1, 3, 4, 7])))
        for i, obs_id in enumerate(sparse_rich.ids(axis='observation')):
            if obs_id == '0':
                self.assertEqual(sparse_rich._observation_metadata[i],
                                 {'Consensus Lineage': 'Bacteria; '
                                  'Actinobacteria; Actinobacteridae; '
                                  'Propionibacterineae; '
                                  'Propionibacterium'})
            elif obs_id == '1':
                self.assertEqual(sparse_rich._observation_metadata[i],
                                 {'Consensus Lineage': 'Bacteria; Firmicutes; '
                                  'Alicyclobacillaceae; Bacilli; '
                                  'Lactobacillales; Lactobacillales; '
                                  'Streptococcaceae; '
                                  'Streptococcus'})
            elif obs_id == '7':
                self.assertEqual(sparse_rich._observation_metadata[i],
                                 {'Consensus Lineage': 'Bacteria; '
                                  'Actinobacteria; Actinobacteridae; '
                                  'Gordoniaceae; '
                                  'Corynebacteriaceae'})
            elif obs_id in ['3', '4']:
                pass  # got lazy
            else:
                raise RuntimeError('obs_id incorrect?')

        self.assertEqual(sparse_rich._sample_metadata, None)

        for i, obs_id in enumerate(sparse_rich.ids(axis='observation')):
            for j, sample_id in enumerate(sparse_rich.ids()):
                if obs_id == '1' and sample_id == 'Key':
                    # should test some abundance data
                    self.assertEqual(True, True)

    def test_from_tsv_dense(self):
        tab1_fh = StringIO(otu_table1)
        sparse_rich = Table.from_tsv(tab1_fh.readlines(), None, None,
                                     OBS_META_TYPES['naive'])
        self.assertEqual(sorted(sparse_rich.ids()),
                         sorted(['Fing', 'Key', 'NA']))
        self.assertEqual(sorted(sparse_rich.ids(axis='observation')),
                         list(map(str, [0, 1, 3, 4, 7])))
        for i, obs_id in enumerate(sparse_rich.ids(axis='observation')):
            if obs_id == '0':
                self.assertEqual(sparse_rich._observation_metadata[i],
                                 {'Consensus Lineage': 'Bacteria; '
                                  'Actinobacteria; Actinobacteridae; '
                                  'Propionibacterineae; '
                                  'Propionibacterium'})
            elif obs_id == '1':
                self.assertEqual(sparse_rich._observation_metadata[i],
                                 {'Consensus Lineage': 'Bacteria; Firmicutes; '
                                  'Alicyclobacillaceae; Bacilli; '
                                  'Lactobacillales; Lactobacillales; '
                                  'Streptococcaceae; '
                                  'Streptococcus'})
            elif obs_id == '7':
                self.assertEqual(sparse_rich._observation_metadata[i],
                                 {'Consensus Lineage': 'Bacteria; '
                                  'Actinobacteria; Actinobacteridae; '
                                  'Gordoniaceae; '
                                  'Corynebacteriaceae'})
            elif obs_id in ['3', '4']:
                pass  # got lazy
            else:
                raise RuntimeError('obs_id incorrect?')

        self.assertEqual(sparse_rich._sample_metadata, None)

        for i, obs_id in enumerate(sparse_rich.ids(axis='observation')):
            for j, sample_id in enumerate(sparse_rich.ids()):
                if obs_id == '1' and sample_id == 'Key':
                    self.assertEqual(True, True)
                    # should test some abundance data

    def test_to_tsv(self):
        """Print out self in a delimited form"""
        exp = '\n'.join(
            ["# Constructed from biom file",
             "#OTU ID\ta\tb",
             "1\t5.0\t6.0",
             "2\t7.0\t8.0"])
        obs = self.st1.delimited_self()
        self.assertEqual(obs, exp)

        # Test observation_column_name.
        exp = '\n'.join(
            ["# Constructed from biom file",
             "Taxon\ta\tb",
             "1\t5.0\t6.0",
             "2\t7.0\t8.0"])
        obs = self.st1.to_tsv(observation_column_name='Taxon')
        self.assertEqual(obs, exp)

    def test_group_metadata_sample(self):
        """Returns the sample group metadata"""
        self.assertEqual(self.st_group_rich.group_metadata(),
                         {'category': ('newick', '(1:0.3,2:0.4);')})
        self.assertEqual(self.st_rich.group_metadata(), None)

    def test_group_metadata_observation(self):
        """Returns the observation group metadata"""
        self.assertEqual(self.st_group_rich.group_metadata(axis='observation'),
                         {'tree': ('newick', '(a:0.3,b:0.4);')})
        self.assertEqual(self.st_rich.group_metadata(axis='observation'), None)

    def test_group_metadata_error(self):
        """Handles invalid input"""
        with self.assertRaises(UnknownAxisError):
            self.simple_derived.group_metadata('bro-axis')

    def test_metadata_invalid_input(self):
        """Correctly handles invalid input."""
        with self.assertRaises(UnknownAxisError):
            self.simple_derived.metadata(1, 'brofist')

    def test_metadata_sample_id(self):
        """returns the sample metadata for a given id"""
        self.assertEqual({'barcode': 'aatt'},
                         self.st_rich.metadata('a'))
        self.assertEqual({'barcode': 'ttgg'},
                         self.st_rich.metadata('b'))

        with self.assertRaises(UnknownIDError):
            self.st_rich.metadata(3, 'sample')

    def test_metadata_sample(self):
        """Return the sample metadata"""
        obs = self.st_rich.metadata()
        exp = [{'barcode': 'aatt'}, {'barcode': 'ttgg'}]
        for o, e in zip(obs, exp):
            self.assertDictEqual(o, e)

    def test_metadata_observation_id(self):
        """returns the observation metadata for a given id"""
        self.assertEqual({'taxonomy': ['k__a', 'p__b']},
                         self.st_rich.metadata('1', 'observation'))
        self.assertEqual({'taxonomy': ['k__a', 'p__c']},
                         self.st_rich.metadata('2', 'observation'))

        with self.assertRaises(UnknownIDError):
            self.simple_derived.metadata('3', 'observation')

    def test_metadata_observation(self):
        """returns the observation metadata"""
        obs = self.st_rich.metadata(axis='observation')
        exp = [{'taxonomy': ['k__a', 'p__b']}, {'taxonomy': ['k__a', 'p__c']}]
        for o, e in zip(obs, exp):
            self.assertDictEqual(o, e)

    def test_index_invalid_input(self):
        """Correctly handles invalid input."""
        with self.assertRaises(UnknownAxisError):
            self.simple_derived.index(1, 'brofist')

    def test_index_sample(self):
        """returns the sample index"""
        self.assertEqual(0, self.simple_derived.index(1, 'sample'))
        self.assertEqual(1, self.simple_derived.index(2, 'sample'))

        with self.assertRaises(UnknownIDError):
            self.simple_derived.index(3, 'sample')

    def test_index_observation(self):
        """returns the observation index"""
        self.assertEqual(0, self.simple_derived.index(3, 'observation'))
        self.assertEqual(1, self.simple_derived.index(4, 'observation'))

        with self.assertRaises(UnknownIDError):
            self.simple_derived.index(5, 'observation')

    def test_index_ids(self):
        """Index all the ids!!!"""
        exp_samp = {1: 0, 2: 1}
        exp_obs = {3: 0, 4: 1}
        self.assertEqual(self.simple_derived._sample_index, exp_samp)
        self.assertEqual(self.simple_derived._obs_index, exp_obs)

    def test_sample_exists(self):
        """Verify samples exist!"""
        self.assertTrue(self.simple_derived.exists(1))
        self.assertTrue(self.simple_derived.exists(2))
        self.assertFalse(self.simple_derived.exists(3))

    def test_observation_exists(self):
        """Verify observation exist!"""
        self.assertTrue(self.simple_derived.exists(3, axis="observation"))
        self.assertTrue(self.simple_derived.exists(4, axis="observation"))
        self.assertFalse(self.simple_derived.exists(2, axis="observation"))

    def test_exists_invalid_axis(self):
        """Verify ValueError raised!"""
        with self.assertRaises(UnknownAxisError):
            self.simple_derived.exists(3, axis="fooz")

    def test_union_id_order(self):
        """Combine unique ids, union"""
        a = [1, 2, 3, 4]
        b = [3, 4, 5, 6, 0, 'a']
        exp = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 0: 6, 'a': 7}
        obs = self.st1._union_id_order(a, b)
        self.assertEqual(obs, exp)

    def test_intersect_id_order(self):
        """Combine ids, intersection"""
        a = [1, 2, 3, 4]
        b = [3, 4, 5, 6, 0, 'a']
        exp = {3: 0, 4: 1}
        obs = self.st1._intersect_id_order(a, b)
        self.assertEqual(obs, exp)

    def test_verify_metadata(self):
        """Make sure the metadata is sane (including obs/sample ids)"""
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 0}, {'e': 0}, {'f': 0}, {'g': 0}]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        Table(d, obs_ids, samp_ids, obs_md, samp_md)
        # test is that no exception is raised

        obs_ids = [1, 2]
        self.assertRaises(TableException, Table, d, samp_ids, obs_ids, samp_md,
                          obs_md)

        obs_ids = [1, 2, 3]
        samp_ids = [4, 5, 6]
        self.assertRaises(TableException, Table, d, samp_ids, obs_ids, samp_md,
                          obs_md)

        samp_ids = [4, 5, 6, 7]
        obs_md = ['a', 'b']
        self.assertRaises(TableException, Table, d, samp_ids, obs_ids, samp_md,
                          obs_md)

        obs_md = ['a', 'b', 'c']
        samp_md = ['d', 'e', 'f']
        self.assertRaises(TableException, Table, d, samp_ids, obs_ids, samp_md,
                          obs_md)

        obs_md = None
        samp_md = None

        # test is that no exception is raised
        Table(d, obs_ids, samp_ids, obs_md, samp_md)

        # do not allow duplicate ids
        obs_ids = [1, 1, 3]
        samp_ids = [4, 5, 6]
        self.assertRaises(TableException, Table, d, samp_ids, obs_ids, samp_md,
                          obs_md)

        obs_ids = [1, 2, 3]
        samp_ids = [4, 4, 6]
        self.assertRaises(TableException, Table, d, samp_ids, obs_ids, samp_md,
                          obs_md)

    def test_cast_metadata(self):
        """Cast metadata objects to defaultdict to support default values"""
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 1}, {'b': 2}, {'c': 3}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 1}, None, {'f': 3}, {'g': 4}]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, obs_md, samp_md)

        self.assertEqual(t._sample_metadata[0]['non existent key'], None)
        self.assertEqual(t._sample_metadata[1]['non existent key'], None)
        self.assertEqual(t._sample_metadata[2]['non existent key'], None)
        self.assertEqual(t._sample_metadata[3]['non existent key'], None)
        self.assertEqual(t._observation_metadata[0]['non existent key'], None)
        self.assertEqual(t._observation_metadata[1]['non existent key'], None)
        self.assertEqual(t._observation_metadata[2]['non existent key'], None)

    def test_add_group_metadata(self):
        """add group metadata works correctly"""
        obs_ids = [1, 2, 3]
        samp_ids = [4, 5, 6, 7]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        obs_g_md = {'tree': ('newick', '(1:0.2,(2:0.3,3:0.4):0.5);')}
        sample_g_md = {'graph': ('edge_list', '(4,5), (4,6), (5,7), (6,7)')}
        t = Table(d, obs_ids, samp_ids, observation_group_metadata=None,
                  sample_group_metadata=sample_g_md)
        t.add_group_metadata(obs_g_md, axis='observation')
        self.assertEqual(t.group_metadata(axis='observation'),
                         {'tree': ('newick', '(1:0.2,(2:0.3,3:0.4):0.5);')})

    def test_add_group_metadata_w_existing_metadata(self):
        """add group metadata works with existing metadata"""
        obs_ids = [1, 2, 3]
        samp_ids = [4, 5, 6, 7]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        obs_g_md = {'tree': ('newick', '(1:0.2,(2:0.3,3:0.4):0.5);')}
        sample_g_md = {'graph': ('edge_list', '(4,5), (4,6), (5,7), (6,7)')}
        t = Table(d, obs_ids, samp_ids, observation_group_metadata=obs_g_md,
                  sample_group_metadata=sample_g_md)
        new_sample_md = {
            'tree': ('newick', '((4:0.1,5:0.1):0.2,(6:0.1,7:0.1):0.2):0.3;')
            }
        t.add_group_metadata(new_sample_md)
        self.assertEqual(
            t.group_metadata(),
            {'graph': ('edge_list', '(4,5), (4,6), (5,7), (6,7)'),
             'tree': ('newick', '((4:0.1,5:0.1):0.2,(6:0.1,7:0.1):0.2):0.3;')})

    def test_to_dataframe(self):
        mat = csr_matrix(np.array([[0.0, 1.0, 2.0],
                                   [3.0, 4.0, 5.0]]))
        exp = pd.DataFrame.sparse.from_spmatrix(mat,
                                                index=['O1', 'O2'],
                                                columns=['S1', 'S2', 'S3'])
        obs = example_table.to_dataframe()

        # assert frame equal between sparse and dense frames wasn't working
        # as expected
        npt.assert_equal(obs.values, exp.values)
        self.assertTrue(all(obs.index == exp.index))
        self.assertTrue(all(obs.columns == exp.columns))

    def test_to_dataframe_is_sparse(self):
        df = example_table.to_dataframe()
        density = (float(example_table.matrix_data.getnnz()) /
                   np.prod(example_table.shape))
        df_density = (df.values > 0).sum().sum() / np.prod(df.shape)
        assert np.allclose(df_density, density)

    def test_to_dataframe_dense(self):
        exp = pd.DataFrame(np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]),
                           index=['O1', 'O2'],
                           columns=['S1', 'S2', 'S3'])
        obs = example_table.to_dataframe(dense=True)
        pdt.assert_frame_equal(obs, exp)

    @pytest.mark.skipif(not HAVE_ANNDATA, reason="anndata not installed")
    def test_to_anndata_dense(self):
        exp = example_table.to_dataframe(dense=True)
        adata = example_table.to_anndata(dense=True, dtype='float64')
        pdt.assert_frame_equal(adata.transpose().to_df(), exp)

    @pytest.mark.skipif(not HAVE_ANNDATA, reason="anndata not installed")
    def test_to_anndata_sparse(self):
        adata = example_table.to_anndata(dense=False)
        mat = example_table.matrix_data.toarray()
        np.testing.assert_array_equal(adata.transpose().X.toarray(), mat)

    @pytest.mark.skipif(not HAVE_ANNDATA, reason="anndata not installed")
    def test_to_anndata_metadata(self):
        adata = example_table.to_anndata()

        obs_samp = example_table.metadata_to_dataframe(axis='sample')
        obs_obs = example_table.metadata_to_dataframe(axis='observation')

        pdt.assert_frame_equal(adata.obs, obs_samp)
        pdt.assert_frame_equal(adata.var, obs_obs)

    def test_metadata_to_dataframe(self):
        exp_samp = pd.DataFrame(['A', 'B', 'A'], index=['S1', 'S2', 'S3'],
                                columns=['environment'])
        exp_obs = pd.DataFrame([['Bacteria', 'Firmicutes'],
                                ['Bacteria', 'Bacteroidetes']],
                               index=['O1', 'O2'],
                               columns=['taxonomy_0', 'taxonomy_1'])
        obs_samp = example_table.metadata_to_dataframe(axis='sample')
        obs_obs = example_table.metadata_to_dataframe(axis='observation')
        pdt.assert_frame_equal(obs_samp, exp_samp)
        pdt.assert_frame_equal(obs_obs, exp_obs)

    def test_metadata_to_dataframe_uneven_list_metadata(self):
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    [{'taxonomy': ['k__foo', 'p__bar']},
                     {'taxonomy': ['k__foo']}])
        exp_obs = pd.DataFrame([['k__foo', 'p__bar'],
                                ['k__foo', None]],
                               index=['a', 'b'],
                               columns=['taxonomy_0', 'taxonomy_1'])
        obs_obs = tab.metadata_to_dataframe(axis='observation')
        pdt.assert_frame_equal(obs_obs, exp_obs)

    def test_metadata_to_dataframe_badaxis(self):
        with self.assertRaises(UnknownAxisError):
            example_table.metadata_to_dataframe(axis='foo')

    def test_metadata_to_dataframe_nomd(self):
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    [{'foo': 1}, {'foo': 2}])
        with self.assertRaises(KeyError):
            tab.metadata_to_dataframe('sample')

    def test_del_metadata_full(self):
        obs_ids = [1, 2, 3]
        obs_md = [{'taxonomy': ['A', 'B'], 'other': 'h1'},
                  {'taxonomy': ['B', 'C'], 'other': 'h2'},
                  {'taxonomy': ['E', 'D', 'F'], 'other': 'h3'}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'foo': 1}, {'foo': 2}, {'foo': 3}, {'foo': 4}]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, observation_metadata=obs_md,
                  sample_metadata=samp_md)
        exp = Table(d, obs_ids, samp_ids)
        t.del_metadata(axis='whole')
        self.assertEqual(t, exp)

    def test_del_metadata_partial(self):
        obs_ids = [1, 2, 3]
        obs_md = [{'taxonomy': ['A', 'B'], 'other': 'h1'},
                  {'taxonomy': ['B', 'C'], 'other': 'h2'},
                  {'taxonomy': ['E', 'D', 'F'], 'other': 'h3'}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 0}, {'e': 0}, {'f': 0}, {'g': 0}]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, observation_metadata=obs_md,
                  sample_metadata=samp_md)

        exp_md = [md.copy() for md in obs_md]
        del exp_md[0]['taxonomy']
        del exp_md[1]['taxonomy']
        del exp_md[2]['taxonomy']

        exp = Table(d, obs_ids, samp_ids, observation_metadata=exp_md,
                    sample_metadata=samp_md)
        t.del_metadata(keys=['taxonomy'], axis='observation')
        self.assertEqual(t, exp)

    def test_del_metadata_nomd(self):
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'])
        exp = tab.copy()
        tab.del_metadata(axis='whole')
        self.assertEqual(tab, exp)

    def test_del_metadata_badaxis(self):
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'])
        with self.assertRaises(UnknownAxisError):
            tab.del_metadata(axis='foo')

    def test_del_metadata_defaults(self):
        tab = example_table.copy()
        tab.del_metadata()
        exp = Table(example_table.matrix_data,
                    example_table.ids(axis='observation'),
                    example_table.ids())
        self.assertEqual(tab, exp)

    def test_del_metadata_idempotent(self):
        ex = example_table.copy()
        ex.del_metadata()
        ex_no_md = ex.copy()
        ex.del_metadata()
        self.assertEqual(ex, ex_no_md)

    def test_del_metadata_empty_list(self):
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    observation_metadata=[{'foo': 1, 'bar': 2, 'baz': 3},
                                          {'foo': 4, 'bar': 5, 'baz': 6}])
        exp = tab.copy()
        tab.del_metadata(keys=[])
        self.assertEqual(tab, exp)

    def test_del_metadata_multiple_keys(self):
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    observation_metadata=[{'foo': 1, 'bar': 2, 'baz': 3},
                                          {'foo': 4, 'bar': 5, 'baz': 6}])
        exp = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    observation_metadata=[{'baz': 3},
                                          {'baz': 6}])
        tab.del_metadata(keys=['foo', 'bar'])
        self.assertEqual(tab, exp)

    def test_del_metadata_jagged(self):
        # this situation should never happen but technically can
        tab = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    observation_metadata=[{'foo': 1}, {'bar': 2}])
        tab.del_metadata(axis='observation', keys=['foo'])
        exp = Table(np.array([[1, 2], [3, 4]]), ['a', 'b'], ['c', 'd'],
                    observation_metadata=[{}, {'bar': 2}])
        self.assertEqual(tab, exp)

    def test_del_metadata_keys_none_sample(self):
        tab = example_table.copy()
        tab.del_metadata(axis='sample')
        self.assertEqual(tab.metadata(), None)

    def test_del_metadata_keys_none_observation(self):
        tab = example_table.copy()
        tab.del_metadata(axis='observation')
        self.assertEqual(tab.metadata(axis='observation'), None)

    def test_del_metadata_keys_none_whole(self):
        tab = example_table.copy()
        tab.del_metadata(axis='whole')
        self.assertEqual(tab.metadata(), None)
        self.assertEqual(tab.metadata(axis='observation'), None)

    def test_all_keys_dropped(self):
        tab = example_table.copy()
        tab.del_metadata(keys=['taxonomy', 'environment'], axis='whole')
        self.assertEqual(tab, Table(tab.matrix_data,
                                    tab.ids(axis='observation'),
                                    tab.ids()))

    def test_add_metadata_two_entries(self):
        """ add_metadata functions with more than one md entry """
        obs_ids = [1, 2, 3]
        obs_md = {1: {'taxonomy': ['A', 'B'], 'other': 'h1'},
                  2: {'taxonomy': ['B', 'C'], 'other': 'h2'},
                  3: {'taxonomy': ['E', 'D', 'F'], 'other': 'h3'}}
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 0}, {'e': 0}, {'f': 0}, {'g': 0}]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, observation_metadata=None,
                  sample_metadata=samp_md)
        t.add_metadata(obs_md, axis='observation')
        self.assertEqual(t._observation_metadata[0]['taxonomy'], ['A', 'B'])
        self.assertEqual(t._observation_metadata[1]['taxonomy'], ['B', 'C'])
        self.assertEqual(t._observation_metadata[2]['taxonomy'],
                         ['E', 'D', 'F'])
        self.assertEqual(t._observation_metadata[0]['other'], 'h1')
        self.assertEqual(t._observation_metadata[1]['other'], 'h2')
        self.assertEqual(t._observation_metadata[2]['other'], 'h3')

        samp_md = {4: {'x': 'y', 'foo': 'bar'}, 5: {'x': 'z'}}
        t.add_metadata(samp_md, axis='sample')
        self.assertEqual(t._sample_metadata[0]['x'], 'y')
        self.assertEqual(t._sample_metadata[0]['foo'], 'bar')
        self.assertEqual(t._sample_metadata[1]['x'], 'z')

    def test_add_metadata_one_w_existing_metadata(self):
        """ add_sample_metadata functions with existing metadata """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'Treatment': 'Control'},
                   {'Treatment': 'Fasting'},
                   {'Treatment': 'Fasting'},
                   {'Treatment': 'Control'}]
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, observation_metadata=obs_md,
                  sample_metadata=samp_md)
        self.assertEqual(t._sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t._sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[3]['Treatment'], 'Control')

        samp_md = {4: {'barcode': 'TTTT'},
                   6: {'barcode': 'AAAA'},
                   5: {'barcode': 'GGGG'},
                   7: {'barcode': 'CCCC'},
                   10: {'ignore': 'me'}}
        t.add_metadata(samp_md, 'sample')
        self.assertEqual(t._sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t._sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[3]['Treatment'], 'Control')
        self.assertEqual(t._sample_metadata[0]['barcode'], 'TTTT')
        self.assertEqual(t._sample_metadata[1]['barcode'], 'GGGG')
        self.assertEqual(t._sample_metadata[2]['barcode'], 'AAAA')
        self.assertEqual(t._sample_metadata[3]['barcode'], 'CCCC')

        obs_md = {1: {'foo': 'bar'}}
        t.add_metadata(obs_md, axis='observation')
        self.assertEqual(t._observation_metadata[0]['foo'], 'bar')
        self.assertEqual(t._observation_metadata[1]['foo'], None)
        self.assertEqual(t._observation_metadata[2]['foo'], None)

    def test_add_metadata_one_entry(self):
        """ add_sample_metadata functions with single md entry """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = {4: {'Treatment': 'Control'},
                   5: {'Treatment': 'Fasting'},
                   6: {'Treatment': 'Fasting'},
                   7: {'Treatment': 'Control'}}
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, obs_md=obs_md, samp_md=None)
        t.add_metadata(samp_md, axis='sample')
        self.assertEqual(t._sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t._sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[3]['Treatment'], 'Control')

    def test_add_sample_metadata_two_entries(self):
        """ add_sample_metadata functions with more than one md entry """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = {4: {'Treatment': 'Control', 'D': ['A', 'A']},
                   5: {'Treatment': 'Fasting', 'D': ['A', 'B']},
                   6: {'Treatment': 'Fasting', 'D': ['A', 'C']},
                   7: {'Treatment': 'Control', 'D': ['A', 'D']}}
        d = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, obs_ids, samp_ids, obs_md=obs_md, samp_md=None)
        t.add_metadata(samp_md, axis='sample')
        self.assertEqual(t._sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t._sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t._sample_metadata[3]['Treatment'], 'Control')
        self.assertEqual(t._sample_metadata[0]['D'], ['A', 'A'])
        self.assertEqual(t._sample_metadata[1]['D'], ['A', 'B'])
        self.assertEqual(t._sample_metadata[2]['D'], ['A', 'C'])
        self.assertEqual(t._sample_metadata[3]['D'], ['A', 'D'])

    def test_align_to_dataframe_samples(self):
        table = Table(np.array([[0, 0, 1, 1],
                                [2, 2, 4, 4],
                                [5, 5, 3, 3],
                                [0, 0, 0, 1]]).T,
                      ['o1', 'o2', 'o3', 'o4'],
                      ['s1', 's2', 's3', 's4'])
        metadata = pd.DataFrame([['a', 'control'],
                                 ['c', 'diseased'],
                                 ['b', 'control']],
                                index=['s1', 's3', 's2'],
                                columns=['Barcode', 'Treatment'])
        exp_table = Table(np.array([[0, 0, 1, 1],
                                    [2, 2, 4, 4],
                                    [5, 5, 3, 3]]).T,
                          ['o1', 'o2', 'o3', 'o4'],
                          ['s1', 's2', 's3'])
        exp_metadata = pd.DataFrame([['a', 'control'],
                                     ['b', 'control'],
                                     ['c', 'diseased']],
                                    index=['s1', 's2', 's3'],
                                    columns=['Barcode', 'Treatment'])
        res_table, res_metadata = table.align_to_dataframe(metadata)
        pdt.assert_frame_equal(exp_metadata, res_metadata)
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')

    def test_align_to_dataframe_observations(self):
        table = Table(np.array([[0, 0, 1, 1],
                                [2, 2, 4, 4],
                                [5, 5, 3, 3],
                                [0, 0, 0, 1]]),
                      ['o1', 'o2', 'o3', 'o4'],
                      ['s1', 's2', 's3', 's4'])
        metadata = pd.DataFrame([['a', 'Firmicutes'],
                                 ['c', 'Proteobacteria'],
                                 ['b', 'Firmicutes']],
                                index=['o1', 'o3', 'o2'],
                                columns=['Barcode', 'Treatment'])
        exp_table = Table(np.array([[0, 0, 1, 1],
                                    [2, 2, 4, 4],
                                    [5, 5, 3, 3]]),
                          ['o1', 'o2', 'o3'],
                          ['s1', 's2', 's3', 's4'])
        exp_metadata = pd.DataFrame([['a', 'Firmicutes'],
                                     ['b', 'Firmicutes'],
                                     ['c', 'Proteobacteria']],
                                    index=['o1', 'o2', 'o3'],
                                    columns=['Barcode', 'Treatment'])
        res_table, res_metadata = table.align_to_dataframe(
            metadata, axis='observation')
        pdt.assert_frame_equal(exp_metadata, res_metadata)
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')

    def test_align_to_dataframe_samples_remove_empty(self):
        table = Table(np.array([[0, 0, 1, 0],
                                [2, 2, 4, 0],
                                [5, 5, 3, 0],
                                [0, 0, 0, 1]]).T,
                      ['o1', 'o2', 'o3', 'o4'],
                      ['s1', 's2', 's3', 's4'])
        metadata = pd.DataFrame([['a', 'control'],
                                 ['c', 'diseased'],
                                 ['b', 'control']],
                                index=['s1', 's3', 's2'],
                                columns=['Barcode', 'Treatment'])
        exp_table = Table(np.array([[0, 0, 1],
                                    [2, 2, 4],
                                    [5, 5, 3]]).T,
                          ['o1', 'o2', 'o3'],
                          ['s1', 's2', 's3'])
        exp_metadata = pd.DataFrame([['a', 'control'],
                                     ['b', 'control'],
                                     ['c', 'diseased']],
                                    index=['s1', 's2', 's3'],
                                    columns=['Barcode', 'Treatment'])
        res_table, res_metadata = table.align_to_dataframe(metadata)
        pdt.assert_frame_equal(exp_metadata, res_metadata)
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')

    def test_align_to_dataframe_samples_empty(self):
        table = Table(np.array([[0, 0, 1, 0],
                                [2, 2, 4, 0],
                                [5, 5, 3, 0],
                                [0, 0, 0, 1]]).T,
                      ['o1', 'o2', 'o3', 'o4'],
                      ['s1', 's2', 's3', 's4'])
        metadata = pd.DataFrame([['a', 'control'],
                                 ['c', 'diseased'],
                                 ['b', 'control']],
                                index=['s1', 's3', 's2'],
                                columns=['Barcode', 'Treatment'])
        exp_table = Table(np.array([[0, 0, 1],
                                    [2, 2, 4],
                                    [5, 5, 3]]).T,
                          ['o1', 'o2', 'o3'],
                          ['s1', 's2', 's3'])
        exp_metadata = pd.DataFrame([['a', 'control'],
                                     ['b', 'control'],
                                     ['c', 'diseased']],
                                    index=['s1', 's2', 's3'],
                                    columns=['Barcode', 'Treatment'])
        res_table, res_metadata = table.align_to_dataframe(metadata)
        pdt.assert_frame_equal(exp_metadata, res_metadata)
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')

    def test_align_to_dataframe_samples_no_common_ids(self):
        table = Table(np.array([[0, 0, 1, 0],
                                [2, 2, 4, 0],
                                [5, 5, 3, 0],
                                [0, 0, 0, 1]]),
                      ['s1', 's2', 's3', 's4'],
                      ['o1', 'o2', 'o3', 'o4'])
        metadata = pd.DataFrame([['a', 'control'],
                                 ['c', 'diseased'],
                                 ['b', 'control']],
                                index=['s1', 's3', 's2'],
                                columns=['Barcode', 'Treatment'])
        with self.assertRaises(TableException):
            table.align_to_dataframe(metadata)

    @pytest.mark.skipif(not HAVE_SKBIO, reason="skbio not installed")
    def test_align_tree_issue_948(self):
        table = Table(np.array([[0, 0, 0, 0],
                                [2, 3, 4, 4],
                                [5, 5, 3, 3],
                                [0, 0, 0, 1]]).T,
                      ['a', 'b', 'c', 'd'],
                      ['s1', 's2', 's3', 's4'])
        tree = skbio.TreeNode.read(["(a,b,c,d)r;"])
        exp_tree = tree
        exp_table = table.copy()
        res_table, res_tree = table.align_tree(tree)
        self.assertEqual(res_table, exp_table)
        self.assertEqual(str(exp_tree), str(res_tree))

    @pytest.mark.skipif(not HAVE_SKBIO, reason="skbio not installed")
    def test_align_tree_intersect_tips(self):
        # there are less tree tips than observations
        table = Table(np.array([[0, 0, 1, 1],
                                [2, 3, 4, 4],
                                [5, 5, 3, 3],
                                [0, 0, 0, 1]]).T,
                      ['a', 'b', 'c', 'd'],
                      ['s1', 's2', 's3', 's4'])
        tree = skbio.TreeNode.read(["((a,b)f,d)r;"])
        exp_table = Table(np.array([[0, 0, 1],
                                    [2, 3, 4],
                                    [5, 5, 3],
                                    [0, 0, 1]]).T,
                          ['a', 'b', 'd'],
                          ['s1', 's2', 's3', 's4'])
        exp_tree = tree
        res_table, res_tree = table.align_tree(tree)
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')
        self.assertEqual(str(exp_tree), str(res_tree))

    @pytest.mark.skipif(not HAVE_SKBIO, reason="skbio not installed")
    def test_align_tree_intersect_obs(self):
        # table has less observations than tree tips
        table = Table(np.array([[0, 0, 1],
                                [2, 3, 4],
                                [5, 5, 3],
                                [0, 0, 1]]).T,
                      ['a', 'b', 'd'],
                      ['s1', 's2', 's3', 's4'])
        tree = skbio.TreeNode.read(["(((a,b)f, c),d)r;"])
        exp_table = Table(np.array([[1, 0, 0],
                                    [4, 2, 3],
                                    [3, 5, 5],
                                    [1, 0, 0]]).T,
                          ['d', 'a', 'b'],
                          ['s1', 's2', 's3', 's4'])
        exp_tree = skbio.TreeNode.read(["(d,(a,b)f)r;"])
        res_table, res_tree = table.align_tree(tree)
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')
        self.assertEqual(exp_tree.compare_rfd(res_tree), 0)

    @pytest.mark.skipif(not HAVE_SKBIO, reason="skbio not installed")
    def test_align_tree_sample(self):
        # table has less observations than tree tips
        table = Table(np.array([[0, 0, 1],
                                [2, 3, 4],
                                [5, 5, 3],
                                [0, 0, 1]]),
                      ['o1', 'o2', 'o3', 'o4'],
                      ['s1', 's2', 's4'])
        tree = skbio.TreeNode.read(["(((s1,s2)F, s3),s4)R;"])
        exp_table = Table(np.array([[1, 0, 0],
                                    [4, 2, 3],
                                    [3, 5, 5],
                                    [1, 0, 0]]),
                          ['o1', 'o2', 'o3', 'o4'],
                          ['s4', 's1', 's2'])
        exp_tree = skbio.TreeNode.read(["(s4,(s1,s2)F)R;"])
        res_table, res_tree = table.align_tree(tree, axis='sample')
        self.assertEqual(res_table.descriptive_equality(exp_table),
                         'Tables appear equal')
        self.assertEqual(exp_tree.compare_rfd(res_tree), 0)

    def test_get_value_by_ids(self):
        """Return the value located in the matrix by the ids"""
        t1 = Table(np.array([[5, 6], [7, 8]]), [3, 4], [1, 2])
        t2 = Table(np.array([[5, 6], [7, 8]]),
                   ['c', 'd'], ['a', 'b'])

        self.assertEqual(5, t1.get_value_by_ids(3, 1))
        self.assertEqual(6, t1.get_value_by_ids(3, 2))
        self.assertEqual(7, t1.get_value_by_ids(4, 1))
        self.assertEqual(8, t1.get_value_by_ids(4, 2))
        self.assertEqual(5, t2.get_value_by_ids('c', 'a'))
        self.assertEqual(6, t2.get_value_by_ids('c', 'b'))
        self.assertEqual(7, t2.get_value_by_ids('d', 'a'))
        self.assertEqual(8, t2.get_value_by_ids('d', 'b'))

        self.assertRaises(UnknownIDError, t1.get_value_by_ids, 'a', 1)
        self.assertRaises(UnknownIDError, t2.get_value_by_ids, 0, 0)

    def test_getitem(self):
        """getitem should work as expeceted"""
        self.assertEqual(self.simple_derived[0, 0], 5)
        self.assertEqual(self.simple_derived[1, 0], 7)
        self.assertEqual(self.simple_derived[0, 1], 6)
        self.assertEqual(self.simple_derived[1, 1], 8)
        self.assertRaises(IndexError, self.simple_derived.__getitem__, [1, 2])

    def test_is_empty(self):
        """returns true if empty"""
        self.assertTrue(Table(np.array([]), [], []).is_empty())
        self.assertFalse(self.simple_derived.is_empty())
        self.assertTrue(Table(np.array([[]]), [], []).is_empty())

    def test_convert_vector_to_dense(self):
        """Properly converts ScipySparseMat vectors to dense numpy repr."""
        input_row = lil_matrix((1, 3))
        input_row[(0, 0)] = 1
        input_row[(0, 2)] = 3
        exp = np.array([1, 0, 3])
        obs = self.row_vec._to_dense(input_row)
        npt.assert_array_equal(obs, exp)

        input_row = lil_matrix((3, 1))
        input_row[(0, 0)] = 1
        input_row[(2, 0)] = 3
        exp = np.array([1, 0, 3])
        obs = self.row_vec._to_dense(input_row)
        npt.assert_array_equal(obs, exp)

        input_row = lil_matrix((1, 1))
        input_row[(0, 0)] = 42
        exp = np.array([42])
        obs = self.single_ele._to_dense(input_row)
        npt.assert_array_equal(obs, exp)

    def test_length(self):
        npt.assert_array_equal(self.null1.length(), 0)
        npt.assert_array_equal(self.null2.length(axis='sample'), 42)
        npt.assert_array_equal(self.null3.length(axis='observation'), 42)
        npt.assert_array_equal(self.mat1.length(), 3)
        npt.assert_array_equal(self.empty.length(axis='observation'), 2)
        npt.assert_array_equal(self.row_vec.length(axis='observation'), 1)
        npt.assert_array_equal(self.row_vec.length(axis='sample'), 3)

        with self.assertRaises(UnknownAxisError):
            self.mat1.length(axis='foo')

    def test_shape(self):
        """What kind of shape are you in?"""
        npt.assert_array_equal(self.null1.shape, (0, 0))
        npt.assert_array_equal(self.null2.shape, (0, 42))
        npt.assert_array_equal(self.null3.shape, (42, 0))
        npt.assert_array_equal(self.mat1.shape, (2, 3))
        npt.assert_array_equal(self.empty.shape, (2, 2))
        npt.assert_array_equal(self.row_vec.shape, (1, 3))
        npt.assert_array_equal(self.col_vec.shape, (3, 1))
        npt.assert_array_equal(self.single_ele.shape, (1, 1))

    def test_dtype(self):
        """What's your type?"""
        for m in self.nulls:
            self.assertEqual(m.dtype, None)

        self.assertEqual(self.empty.dtype, float)
        self.assertEqual(self.row_vec.dtype, float)

    def test_nnz(self):
        """What is your NNZ?"""
        for m in self.nulls:
            self.assertEqual(m.nnz, 0)

        self.assertEqual(self.empty.nnz, 0)
        self.assertEqual(self.single_ele.nnz, 1)
        self.assertEqual(self.mat1.nnz, 4)
        self.assertEqual(self.explicit_zeros.nnz, 4)

    def test_nnz_issue_727(self):
        wrn = "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient."  # noqa
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=wrn)
            tab = Table(np.array([[0, 1], [0, 0]]), ['a', 'b'], ['1', '2'])
            self.assertEqual(tab.nnz, 1)
            tab._data[0, 0] = 0
            self.assertEqual(tab.nnz, 1)

    def test_get_row(self):
        """Test grabbing a row from the matrix."""
        # note that we only have to test the first two elements don't have that
        # row according to the underlying scipy sparse matrix
        for i in range(0, 2):
            with self.assertRaises(IndexError):
                self.nulls[i]._get_row(0)

        exp = lil_matrix((1, 3))
        exp[(0, 0)] = 1
        exp[(0, 2)] = 2

        obs = self.mat1._get_row(0)
        self.assertEqual((obs != exp).sum(), 0)

    def test_get_col(self):
        """Test grabbing a column from the matrix."""
        # note that we only have to test the first and last element, these
        # don't have that column according to the underlying scipy sparse
        # matrix
        for i in [0, 2]:
            with self.assertRaises(IndexError):
                self.nulls[i]._get_col(0)

        exp = lil_matrix((2, 1))
        exp[(0, 0)] = 1
        exp[(1, 0)] = 3

        obs = self.mat1._get_col(0)
        self.assertEqual((obs != exp).sum(), 0)

    def test_eq(self):
        """Test whether two matrices are equal."""
        # Empty/null cases (i.e., 0x0, 0xn, nx0).
        def ids(X):
            return ['x%d' % e for e in range(0, X)]
        a = Table(np.zeros((0, 0)), [], [])
        b = Table(np.zeros((0, 42), dtype=float), [], ids(42))
        c = Table(np.zeros((42, 0), dtype=float), ids(42), [])
        d = Table(np.zeros((2, 2)), ids(2), ids(2))

        self.assertTrue(self.null1 == a)
        self.assertTrue(self.null2 == b)
        self.assertTrue(self.null3 == c)
        self.assertTrue(self.empty == d)

        mat2 = Table(np.array([[1, 0, 2], [3, 0, 4]]),
                     ['o1', 'o2'], ['s1', 's2', 's3'])
        self.assertTrue(self.mat1 == mat2)

        mat2._data = mat2._data.tolil()
        self.assertNotEqual(self.mat1._data.format, mat2._data.format)
        self.assertEqual(self.mat1, mat2)

        # Equality works in both directions.
        self.assertEqual(mat2, self.mat1)

    def test_ne(self):
        """Test whether two matrices are not equal."""
        # Wrong type.
        self.assertTrue(self.null1 != np.array([]))

        # Wrong shape.
        def ids(X):
            return ['x%d' % e for e in range(0, X)]
        d = Table(np.ones((1, 1)), ids(1), ids(1))
        self.assertTrue(self.null2 != self.null3)
        self.assertTrue(self.empty != d)

        # Wrong dtype.
        d = Table(np.zeros((2, 2)), ids(2), ids(2), type=float)
        self.assertTrue(self.empty != d)

        # Wrong size.
        wrong_size = Table(np.zeros((2, 2)), ids(2), ids(2))
        self.assertTrue(self.empty == wrong_size)
        wrong_size = Table(np.ones((1, 1)), ['c'], ['a'])
        self.assertTrue(self.empty != wrong_size)

        # Wrong size.
        wrong_data = self.mat1.copy()
        self.assertTrue(self.mat1 == wrong_data)
        wrong_data = Table(np.array([[42, 0, 2], [3, 0, 4]]),
                           ['o1', 'o2'], ['s1', 's2', 's3'])
        self.assertTrue(self.mat1 != wrong_data)
        self.assertTrue(wrong_data != self.mat1)

    def test_getitem_2(self):
        """Test getting an element from the matrix."""
        for m in self.nulls:
            with self.assertRaises(IndexError):
                m[0, 0]

        with self.assertRaises(IndexError):
            self.empty[0]

        with self.assertRaises(IndexError):
            self.empty[:, :]

        with self.assertRaises(IndexError):
            self.empty[0:1, 0]

        with self.assertRaises(IndexError):
            self.empty[0, 0:1]

        exp = lil_matrix((2, 1))
        obs = self.empty[:, 0]
        self.assertEqual((obs != exp).sum(), 0)

        # Extracting a column.
        obs = self.mat1[:, 2]
        self.assertEqual((obs != self.mat1._get_col(2)).sum(), 0)

        # Extracting a row.
        obs = self.mat1[1, :]
        self.assertEqual((obs != self.mat1._get_row(1)).sum(), 0)

        # Extracting a single element.
        self.assertEqual(self.empty[1, 1], 0)
        self.assertEqual(self.mat1[1, 2], 4)

        with self.assertRaises(IndexError):
            self.mat1[1, 3]


class SparseTableTests(TestCase):

    def setUp(self):
        self.vals = {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8}
        self.st1 = Table(self.vals,
                         ['1', '2'], ['a', 'b'])
        self.st2 = Table(self.vals,
                         ['1', '2'], ['a', 'b'])
        self.vals3 = {(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}
        self.vals4 = {(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}
        self.st3 = Table(self.vals3, ['2', '3'], ['b', 'c'])
        self.st4 = Table(self.vals4, ['3', '4'], ['c', 'd'])
        self._to_dict_f = lambda x: sorted(x.items())
        self.st_rich = Table(self.vals,
                             ['1', '2'], ['a', 'b'],
                             [{'taxonomy': ['k__a', 'p__b']},
                              {'taxonomy': ['k__a', 'p__c']}],
                             [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])

        self.empty_st = Table([], [], [])

        self.vals5 = {(0, 1): 2, (1, 1): 4}
        self.st5 = Table(self.vals5, ['5', '6'], ['a', 'b'])

        self.vals6 = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
        self.st6 = Table(self.vals6, ['5', '6'], ['a', 'b'])

        self.vals7 = {(0, 0): 5, (0, 1): 7, (1, 0): 8, (1, 1): 0}
        self.st7 = Table(self.vals7, ['5', '6'], ['a', 'b'])

        self.single_sample_st = Table(
            np.array([[2.0], [0.0], [1.0]]),
            ['O1', 'O2', 'O3'], ['S1'])
        self.single_obs_st = Table(np.array([[2.0, 0.0, 1.0]]),
                                   ['01'], ['S1', 'S2', 'S3'])

        self.sparse_table = Table(np.array([[1, 0, 2, 0],
                                            [0, 3, 4, 0],
                                            [0, 5, 0, 0]]),
                                  ['O1', 'O2', 'O3'],
                                  ['S1', 'S2', 'S3', 'S4'])

    def test_sum(self):
        """Test of sum!"""
        self.assertEqual(self.st1.sum('whole'), 26)
        npt.assert_equal(self.st1.sum('sample'), np.array([12, 14]))
        npt.assert_equal(self.st1.sum('observation'), np.array([11, 15]))

        exp = np.array([3.0])
        obs = self.single_sample_st.sum('sample')
        self.assertEqual(obs, exp)

        exp = np.array([3.0])
        obs = self.single_obs_st.sum('observation')
        self.assertEqual(obs, exp)

    def test_reduce(self):
        """Reduce method"""
        def f(x, y):
            return x * 2 + y
        npt.assert_equal(self.st1.reduce(f, 'sample'), np.array([17, 20]))
        npt.assert_equal(self.st1.reduce(f, 'observation'), np.array([16, 22]))

    def test_transpose(self):
        """Should transpose a sparse table"""
        obs = self.st1.transpose()

        npt.assert_equal(obs.ids(), self.st1.ids(axis='observation'))
        npt.assert_equal(obs.ids(axis='observation'), self.st1.ids())
        npt.assert_equal(obs.data('1', 'sample'),
                         self.st1.data('1', 'observation'))
        npt.assert_equal(obs.data('2', 'sample'),
                         self.st1.data('2', 'observation'))
        self.assertEqual(obs.transpose(), self.st1)

        obs = self.st_rich.transpose()

        npt.assert_equal(obs.ids(), self.st_rich.ids(axis='observation'))
        npt.assert_equal(obs.ids(axis='observation'), self.st_rich.ids())
        self.assertEqual(obs._sample_metadata,
                         self.st_rich._observation_metadata)
        self.assertEqual(obs._observation_metadata,
                         self.st_rich._sample_metadata)
        npt.assert_equal(obs.data('1', 'sample'),
                         self.st_rich.data('1', 'observation'))
        npt.assert_equal(obs.data('2', 'sample'),
                         self.st_rich.data('2', 'observation'))
        self.assertEqual(obs.transpose(), self.st_rich)

    def test_update_ids_strict_dtype_bug_issue_957(self):
        t = Table(np.arange(6).reshape(2, 3),
                  ['O1', 'O2'],
                  ['ab', 'cdef', 'ghijkl'])
        exp = Table(np.arange(6).reshape(2, 3),
                    ['O1', 'O2'],
                    ['AB', 'cdef', 'ghijkl'])
        obs = t.update_ids({'ab': 'AB'}, strict=False, inplace=False)
        self.assertEqual(obs, exp)

    def test_update_ids_inplace_bug_892(self):
        t = example_table.copy()
        exp = t.ids().copy()
        with self.assertRaises(TableException):
            t.update_ids({i: 'foo' for i in t.ids()}, inplace=True)
        npt.assert_equal(t.ids(), exp)

    def test_update_ids(self):
        """ids are updated as expected"""
        # update observation ids
        exp = self.st1.copy()
        exp._observation_ids = np.array(['41', '42long'])
        id_map = {'2': '42long', '1': '41'}
        obs = self.st1.update_ids(id_map, axis='observation', inplace=False)
        self.assertEqual(obs, exp)

        # update sample ids
        exp = self.st1.copy()
        exp._sample_ids = np.array(['99', '100'])
        id_map = {'a': '99', 'b': '100'}
        obs = self.st1.update_ids(id_map, axis='sample', inplace=False)
        self.assertEqual(obs, exp)

        # extra ids in id_map are ignored
        exp = self.st1.copy()
        exp._observation_ids = np.array(['41', '42'])
        id_map = {'2': '42', '1': '41', '0': '40'}
        obs = self.st1.update_ids(id_map, axis='observation', inplace=False)
        self.assertEqual(obs, exp)

        # missing ids in id_map when strict=True
        with self.assertRaises(TableException):
            self.st1.update_ids({'b': '100'}, axis='sample', strict=True,
                                inplace=False)

        # missing ids in id_map when strict=False
        exp = self.st1.copy()
        exp._sample_ids = np.array(['a', '100'])
        id_map = {'b': '100'}
        obs = self.st1.update_ids(id_map, axis='sample', strict=False,
                                  inplace=False)
        self.assertEqual(obs, exp)

        # raise an error if update would result in duplicated ids
        with self.assertRaises(TableException):
            self.st1.update_ids({'a': '100', 'b': '100'}, axis='sample',
                                inplace=False)

        # raises an error if a invalid axis is passed
        with self.assertRaises(UnknownAxisError):
            self.st1.update_ids(id_map, axis='foo', inplace=False)

        # when inplace == False, the input object is unchanged
        exp = self.st1.copy()
        exp._observation_ids = np.array(['41', '42'])
        id_map = {'2': '42', '1': '41'}
        obs = self.st1.update_ids(id_map, axis='observation', inplace=False)
        npt.assert_equal(self.st1._observation_ids, np.array(['1', '2']))
        # when inplace == True, the input object is changed
        obs = self.st1.update_ids(id_map, axis='observation', inplace=True)
        npt.assert_equal(self.st1._observation_ids, np.array(['41', '42']))

    def test_update_ids_nochange_bug(self):
        """ids are updated as expected"""
        # update observation ids
        exp = self.st1.copy()
        id_map = {'1': '1', '2': '2'}
        obs = self.st1.update_ids(id_map, axis='observation', inplace=False)
        self.assertEqual(obs, exp)

        # test having one ID remain unchanged
        exp = self.st1.copy()
        exp._observation_ids = np.array(['1', '3'])
        id_map = {'1': '1', '2': '3'}
        obs = self.st1.update_ids(id_map, axis='observation', inplace=False)
        self.assertEqual(obs, exp)

    def test_update_ids_cache_bug(self):
        obs = self.st1.update_ids({'1': 'x', '2': 'y'}, axis='observation',
                                  inplace=False)
        exp_index = {'x': 0, 'y': 1}
        self.assertEqual(obs._obs_index, exp_index)

        obs = self.st1.update_ids({'a': 'x', 'b': 'y'}, inplace=False)
        exp_index = {'x': 0, 'y': 1}
        self.assertEqual(obs._sample_index, exp_index)

    def test_other_spmatrix_type(self):
        ss = scipy.sparse
        for c in [ss.lil_matrix, ss.bsr_matrix, ss.coo_matrix, ss.dia_matrix,
                  ss.dok_matrix, ss.csc_matrix, ss.csr_matrix]:
            mat = c((2, 2))
            t = Table(mat, ['a', 'b'], [1, 2])
            self.assertTrue(isinstance(t.matrix_data,
                                       (csr_matrix, csc_matrix)))

    def test_sort_order(self):
        """sorts tables by arbitrary order"""
        # sort by observations arbitrary order
        vals = {(0, 0): 7, (0, 1): 8, (1, 0): 5, (1, 1): 6}
        exp = Table(vals, ['2', '1'], ['a', 'b'])
        obs = self.st1.sort_order(['2', '1'], axis='observation')
        self.assertEqual(obs, exp)
        # sort by samples arbitrary order
        vals = {(0, 0): 6, (0, 1): 5,
                (1, 0): 8, (1, 1): 7}
        exp = Table(vals, ['1', '2'], ['b', 'a'])
        obs = self.st1.sort_order(['b', 'a'], axis='sample')
        self.assertEqual(obs, exp)
        # raises an error if a invalid axis is passed
        with self.assertRaises(UnknownAxisError):
            self.st1.sort_order(['b', 'a'], axis='foo')

    def test_sort(self):
        """table sorted by a function and provided axis"""
        # sort by samples by a function
        sort_f = sorted
        data_in = nparray_to_sparse(
            np.array([[1, 2, 3, 8], [4, 5, 6, 9], [7, 8, 9, 11]]))
        t = Table(data_in, [2, 1, 3], ['c', 'a', 'b', 'd'])
        exp_data = nparray_to_sparse(
            np.array([[2, 3, 1, 8], [5, 6, 4, 9], [8, 9, 7, 11]]))
        exp = Table(exp_data, [2, 1, 3], ['a', 'b', 'c', 'd'])
        obs = t.sort(sort_f=sort_f)
        self.assertEqual(obs, exp)
        # sort by observation ids by a function
        sort_f = sorted
        data_in = nparray_to_sparse(
            np.array([[1, 2, 3, 8], [4, 5, 6, 9], [7, 8, 9, 11]]), float)
        t = Table(data_in, [2, 1, 3], ['c', 'a', 'b', 'd'])
        exp_data = nparray_to_sparse(
            np.array([[4, 5, 6, 9], [1, 2, 3, 8], [7, 8, 9, 11]]), float)
        exp = Table(exp_data, [1, 2, 3], ['c', 'a', 'b', 'd'])
        obs = t.sort(sort_f=sort_f, axis='observation')
        self.assertEqual(obs, exp)
        # raises an error if a invalid axis is passed
        with self.assertRaises(UnknownAxisError):
            t.sort(axis='foo')

    def test_eq(self):
        """sparse equality"""
        self.assertTrue(self.st1 == self.st2)
        self.st1._observation_ids = np.array(["1", "2", "3"], dtype=object)
        self.assertFalse(self.st1 == self.st2)

        self.st1._observation_ids = self.st2._observation_ids
        self.st1._data = nparray_to_sparse(np.array([[1, 2], [10, 20]]))
        self.assertFalse(self.st1 == self.st2)

    def test_data_equality(self):
        """check equality between tables"""
        self.assertTrue(self.st1._data_equality(self.st2._data))
        self.assertTrue(self.st1._data_equality(self.st1._data))
        self.assertFalse(self.st1._data_equality(self.st3._data))

    def test_nonzero(self):
        """Return a list of nonzero positions"""
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(data, ['1', '2', '3'], ['a', 'b', 'c', 'd'])
        exp = [('1', 'a'), ('1', 'b'), ('1', 'd'), ('2', 'b'), ('2', 'd'),
               ('3', 'a'), ('3', 'b')]
        obs = list(st.nonzero())
        self.assertEqual(obs, exp)

    def test_nonzero_csc_bug(self):
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(data, ['1', '2', '3'], ['a', 'b', 'c', 'd'])
        st._data = st._data.tocsc()
        exp = [('1', 'a'), ('1', 'b'), ('1', 'd'), ('2', 'b'), ('2', 'd'),
               ('3', 'a'), ('3', 'b')]
        obs = list(st.nonzero())
        self.assertEqual(obs, exp)

    def test_nonzero_counts(self):
        """Returns nonzero counts over an axis"""
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(data, ['1', '2', '3'], ['a', 'b', 'c', 'd'])

        exp_samp = np.array([6, 12, 0, 11])
        exp_obs = np.array([14, 15, 0])
        exp_whole = np.array([29])

        obs_samp = st.nonzero_counts('sample', binary=False)
        obs_obs = st.nonzero_counts('observation', binary=False)
        obs_whole = st.nonzero_counts('whole', binary=False)

        npt.assert_equal(obs_samp, exp_samp)
        npt.assert_equal(obs_obs, exp_obs)
        npt.assert_equal(obs_whole, exp_whole)

    def test_nonzero_counts_binary(self):
        """Returns nonzero counts over an axis"""
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(data, ['1', '2', '3'], ['a', 'b', 'c', 'd'])

        exp_samp = np.array([2, 3, 0, 2])
        exp_obs = np.array([3, 2, 2])
        exp_whole = np.array([7])

        obs_samp = st.nonzero_counts('sample', binary=True)
        obs_obs = st.nonzero_counts('observation', binary=True)
        obs_whole = st.nonzero_counts('whole', binary=True)

        npt.assert_equal(obs_samp, exp_samp)
        npt.assert_equal(obs_obs, exp_obs)
        npt.assert_equal(obs_whole, exp_whole)

    def test_fast_merge(self):
        data = {(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16}
        exp = Table(data, ['1', '2'], ['a', 'b'])
        obs = self.st1._fast_merge([self.st1])
        self.assertEqual(obs, exp)

    def test_fast_merge_multiple(self):
        data = {(0, 0): 20, (0, 1): 24, (1, 0): 28, (1, 1): 32}
        exp = Table(data, ['1', '2'], ['a', 'b'])
        obs = self.st1._fast_merge([self.st1, self.st1, self.st1])
        self.assertEqual(obs, exp)

    def test_fast_merge_nonoverlapping(self):
        t2 = self.st1.copy()
        t2.update_ids({'a': 'd'}, inplace=True, strict=False)
        t2.update_ids({'2': '3'}, axis='observation', inplace=True,
                      strict=False)
        exp = Table(np.array([[5, 12, 5],
                              [7, 8, 0],
                              [0, 8, 7]]), ['1', '2', '3'],
                    ['a', 'b', 'd'])
        obs = t2._fast_merge([self.st1])
        self.assertEqual(obs, exp)

    def test_merge(self):
        """Merge two tables"""
        u = 'union'
        i = 'intersection'

        # test 1
        data = {(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16}
        exp = Table(data, ['1', '2'], ['a', 'b'])
        obs = self.st1.merge(self.st1, sample=u, observation=u)
        self.assertEqual(obs, exp)

        # test 2
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (1, 0): 7, (1, 1): 9,
                (1, 2): 2, (2, 0): 0, (2, 1): 3, (2, 2): 4}
        exp = Table(data, ['1', '2', '3'], ['a', 'b', 'c'])
        obs = self.st1.merge(self.st3, sample=u, observation=u)
        self.assertEqual(obs, exp)

        # test 3
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 0,
                (1, 0): 7, (1, 1): 8, (1, 2): 0, (1, 3): 0,
                (2, 0): 0, (2, 1): 0, (2, 2): 1, (2, 3): 2,
                (3, 0): 0, (3, 1): 0, (3, 2): 3, (3, 3): 4}
        exp = Table(data, ['1', '2', '3', '4'], ['a', 'b', 'c', 'd'])
        obs = self.st1.merge(self.st4, sample=u, observation=u)
        self.assertEqual(obs, exp)

        # test 4
        data = {(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16}
        exp = Table(data, ['1', '2'], ['a', 'b'])
        obs = self.st1.merge(self.st1, sample=i, observation=i)
        self.assertEqual(obs, exp)

        # test 5
        exp = Table({(0, 0): 9}, ['2'], ['b'])
        obs = self.st1.merge(self.st3, sample=i, observation=i)
        self.assertEqual(obs, exp)

        # test 6
        self.assertRaises(TableException, self.st1.merge, self.st4, i, i)

        # test 7
        data = {(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16}
        exp = Table(data, ['1', '2'], ['a', 'b'])
        obs = self.st1.merge(self.st1, sample=i, observation=u)
        self.assertEqual(obs, exp)

        # test 8
        data = {(0, 0): 6, (1, 0): 9, (2, 0): 3}
        exp = Table(data, ['1', '2', '3'], ['b'])
        obs = self.st1.merge(self.st3, sample=i, observation=u)
        self.assertEqual(obs, exp)

        # test 9
        self.assertRaises(TableException, self.st1.merge, self.st4, i, u)

        # test 10
        data = {(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16}
        exp = Table(data, ['1', '2'], ['a', 'b'])
        obs = self.st1.merge(self.st1, sample=u, observation=i)
        self.assertEqual(obs, exp)

        # test 11
        data = {(0, 0): 7, (0, 1): 9, (0, 2): 2}
        exp = Table(data, ['2'], ['a', 'b', 'c'])
        obs = self.st1.merge(self.st3, sample=u, observation=i)
        self.assertEqual(obs, exp)

        # test 12
        self.assertRaises(TableException, self.st1.merge, self.st4, u, i)

    def test_data(self):
        """"""
        # Returns observations for a given sample
        exp = np.array([5, 7])
        obs = self.st1.data('a', 'sample')
        npt.assert_equal(obs, exp)
        with self.assertRaises(UnknownIDError):
            self.st1.data('asdasd', 'sample')

        # Returns samples for a given observation
        exp = np.array([5, 6])
        obs = self.st1.data('1', 'observation')
        npt.assert_equal(obs, exp)
        with self.assertRaises(UnknownIDError):
            self.st1.data('asdsad', 'observation')

        # Raises an error with unknown axis
        with self.assertRaises(UnknownAxisError):
            obs = self.st1.data('a', axis='foo')

    def test_data_sparse(self):
        # Returns observations for a given sample
        exp = csc_matrix(np.array([[5], [7]]))
        obs = self.st1.data('a', 'sample', dense=False)
        self.assertEqual((obs != exp).nnz, 0)
        with self.assertRaises(UnknownIDError):
            self.st1.data('asdasd', 'sample')

        # Returns samples for a given observation
        exp = csr_matrix(np.array([5, 6]))
        obs = self.st1.data('1', 'observation', dense=False)
        self.assertEqual((obs != exp).nnz, 0)
        with self.assertRaises(UnknownIDError):
            self.st1.data('asdsad', 'observation')

    def test_delimited_self(self):
        """Print out self in a delimited form"""
        exp = '\n'.join(
            ["# Constructed from biom file",
             "#OTU ID\ta\tb",
             "1\t5.0\t6.0",
             "2\t7.0\t8.0"])
        obs = self.st1.delimited_self()
        self.assertEqual(obs, exp)

        # Test observation_column_name.
        exp = '\n'.join(
            ["# Constructed from biom file",
             "Taxon\ta\tb",
             "1\t5.0\t6.0",
             "2\t7.0\t8.0"])
        obs = self.st1.delimited_self(observation_column_name='Taxon')
        self.assertEqual(obs, exp)

    def test_conv_to_self_type(self):
        """Should convert other to sparse type"""
        exp = lil_matrix((2, 2))
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(1, 0)] = 7
        exp[(1, 1)] = 8
        obs = self.st1._conv_to_self_type(self.vals)
        self.assertEqual((obs != exp).sum(), 0)

        exp = lil_matrix((2, 2))
        exp[(0, 0)] = 5
        exp[(0, 1)] = 7
        exp[(1, 0)] = 6
        exp[(1, 1)] = 8
        obs = self.st1._conv_to_self_type(self.vals, transpose=True)
        self.assertEqual((obs != exp).sum(), 0)

        # passing a single vector
        exp = lil_matrix((1, 3))
        exp[(0, 0)] = 2
        exp[(0, 1)] = 0
        exp[(0, 2)] = 3
        obs = self.st1._conv_to_self_type(np.array([2, 0, 3]))
        self.assertEqual((obs != exp).sum(), 0)

        # passing a list of dicts
        exp = lil_matrix((2, 3))
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(0, 2)] = 7
        exp[(1, 0)] = 8
        exp[(1, 1)] = 9
        exp[(1, 2)] = 10
        obs = self.st1._conv_to_self_type([{(0, 0): 5, (0, 1): 6, (0, 2): 7},
                                           {(1, 0): 8, (1, 1): 9, (1, 2): 10}])
        self.assertEqual((obs != exp).sum(), 0)

    def test_to_dense(self):
        """Should convert a self styled vector to numpy type"""
        input_row = lil_matrix((1, 3))
        input_row[(0, 0)] = 10
        exp = np.array([10.0, 0, 0])
        obs = self.st1._to_dense(input_row)
        npt.assert_equal(obs, exp)

        input_col = lil_matrix((3, 1))
        input_col[(0, 0)] = 12
        exp = np.array([12.0, 0, 0])
        obs = self.st1._to_dense(input_col)
        npt.assert_equal(obs, exp)

        # 1x1
        input_vec = lil_matrix((1, 1))
        input_vec[(0, 0)] = 42
        exp = np.array([42.0])
        obs = self.st1._to_dense(input_vec)
        npt.assert_equal(obs, exp)

    def test_iter_data_dense(self):
        exp = [np.array([5, 7]), np.array([6, 8])]
        obs = list(self.st1.iter_data())
        npt.assert_equal(obs, exp)

    def test_iter_data_sparse(self):
        exp = [csr_matrix(np.array([5, 7])),
               csr_matrix(np.array([6, 8]))]
        obs = list(self.st1.iter_data(dense=False))
        for o, e in zip(obs, exp):
            self.assertTrue((o != e).nnz == 0)

    def test_iter_pairwise_simple(self):
        """Should iterate pairwise over samples"""
        exp = [((np.array([5, 7]), 'a', None), (np.array([5, 7]), 'a', None)),
               ((np.array([5, 7]), 'a', None), (np.array([6, 8]), 'b', None)),
               ((np.array([6, 8]), 'b', None), (np.array([5, 7]), 'a', None)),
               ((np.array([6, 8]), 'b', None), (np.array([6, 8]), 'b', None))]
        obs = list(self.st1.iter_pairwise(dense=True, tri=False, diag=True))
        npt.assert_equal(obs, exp)

    def test_iter_pairwise_tri(self):
        """Should iterate pairwise over samples"""
        exp = [((np.array([5, 7]), 'a', None), (np.array([5, 7]), 'a', None)),
               ((np.array([5, 7]), 'a', None), (np.array([6, 8]), 'b', None)),
               ((np.array([6, 8]), 'b', None), (np.array([6, 8]), 'b', None))]
        obs = list(self.st1.iter_pairwise(dense=True, tri=True, diag=True))
        npt.assert_equal(obs, exp)

    def test_iter_pairwise_tri_diag(self):
        """Should iterate pairwise over samples"""
        exp = [((np.array([5, 7]), 'a', None), (np.array([6, 8]), 'b', None))]
        obs = list(self.st1.iter_pairwise(dense=True, tri=True, diag=False))
        npt.assert_equal(obs, exp)

    def test_iter_pairwise_diag(self):
        """Should iterate pairwise over samples"""
        exp = [((np.array([5, 7]), 'a', None), (np.array([6, 8]), 'b', None)),
               ((np.array([6, 8]), 'b', None), (np.array([5, 7]), 'a', None))]
        obs = list(self.st1.iter_pairwise(dense=True, tri=False, diag=False))
        npt.assert_equal(obs, exp)

    def test_iter(self):
        """Should iterate over samples"""
        exp = [(np.array([5, 7]), 'a', None), (np.array([6, 8]), 'b', None)]
        obs = list(self.st1)
        npt.assert_equal(obs, exp)

    def test_iter_obs(self):
        """Iterate over observations of sparse matrix"""
        r1 = lil_matrix((1, 2))
        r2 = lil_matrix((1, 2))
        r1[(0, 0)] = 5
        r1[(0, 1)] = 6
        r2[(0, 0)] = 7
        r2[(0, 1)] = 8

        exp = [r1.tocsr(), r2.tocsr()]
        obs = list(self.st1._iter_obs())

        for o, e in zip(obs, exp):
            self.assertEqual((o != e).sum(), 0)

    def test_iter_samp(self):
        """Iterate over samples of sparse matrix"""
        c1 = lil_matrix((1, 2))
        c2 = lil_matrix((1, 2))
        c1[(0, 0)] = 5
        c1[(0, 1)] = 7
        c2[(0, 0)] = 6
        c2[(0, 1)] = 8

        exp = [c1.tocsc(), c2.tocsc()]
        obs = list(self.st1._iter_samp())

        for o, e in zip(obs, exp):
            self.assertEqual((o != e).sum(), 0)

    def test_iter_samples(self):
        """Iterates samples"""
        gen = self.st1.iter()
        exp = [(np.array([5, 7]), 'a', None), (np.array([6, 8]), 'b', None)]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter()
        exp = [(np.array([5, 7]), 'a', {'barcode': 'aatt'}),
               (np.array([6, 8]), 'b', {'barcode': 'ttgg'})]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        # [[1,2,3],[1,0,2]] isn't yielding column 2 correctly
        vals = {(0, 0): 5, (0, 1): 6, (1, 1): 8}
        st = Table(vals, ['1', '2'], ['a', 'b'])
        gen = st.iter(axis='sample')
        exp = [(np.array([5, 0]), 'a', None), (np.array([6, 8]), 'b', None)]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_observations(self):
        """Iterates observations"""
        gen = self.st1.iter(axis='observation')
        exp = [(np.array([5, 6]), '1', None), (np.array([7, 8]), '2', None)]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter(axis='observation')
        exp = [(np.array([5, 6]), '1', {'taxonomy': ['k__a', 'p__b']}),
               (np.array([7, 8]), '2', {'taxonomy': ['k__a', 'p__c']})]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_sample_data(self):
        """Iterates data by samples"""
        gen = self.st1.iter_data()
        exp = [np.array([5, 7]), np.array([6, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter_data()
        exp = [np.array([5, 7]), np.array([6, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        # [[1,2,3],[1,0,2]] isn't yielding column 2 correctly
        vals = {(0, 0): 5, (0, 1): 6, (1, 1): 8}
        st = Table(vals, ['1', '2'], ['a', 'b'])
        gen = st.iter_data()
        exp = [np.array([5, 0]), np.array([6, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_sample_data_single_obs(self):
        """Iterates data by samples with a single observation."""
        exp = [np.array([2.0]), np.array([0.0]), np.array([1.0])]
        obs = list(self.single_obs_st.iter_data())
        # We test this way to make sure the observed value is a single element
        # array instead of a numpy scalar.
        for o, e in zip(obs, exp):
            self.assertEqual(o, e)

    def test_iter_observation_data(self):
        """Iterates data by observations"""
        gen = self.st1.iter_data(axis="observation")
        exp = [np.array([5, 6]), np.array([7, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter_data(axis="observation")
        exp = [np.array([5, 6]), np.array([7, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_observation_data_single_sample(self):
        """Iterates data by observations from a single sample."""
        exp = [np.array([2.0]), np.array([0.0]), np.array([1.0])]
        obs = list(self.single_sample_st.iter_data(axis="observation"))
        for o, e in zip(obs, exp):
            self.assertEqual(o, e)

    def test_copy_metadata(self):
        self.st_rich._sample_metadata[0]['foo'] = ['bar']
        copied_table = self.st_rich.copy()
        copied_table._sample_metadata[0]['foo'].append('bar2')
        self.assertNotEqual(copied_table, self.st_rich)
        self.st_rich._observation_metadata[0]['foo'] = ['bar']
        copied_table = self.st_rich.copy()
        copied_table._observation_metadata[0]['foo'].append('bar2')
        self.assertNotEqual(copied_table, self.st_rich)

    def test_copy_ids(self):
        copied_table = self.st_rich.copy()
        self.st_rich._sample_ids[0] = 'X'
        self.assertNotEqual(copied_table, self.st_rich)
        copied_table = self.st_rich.copy()
        self.st_rich._observation_ids[0] = 'X'
        self.assertNotEqual(copied_table, self.st_rich)

    def test_copy_data(self):
        copied_table = self.st_rich.copy()
        self.st_rich._data *= 2
        self.assertNotEqual(copied_table, self.st_rich)

    def test_filter_table_with_zeros(self):
        table = self.sparse_table

        def f_sample(vals, id_, md):
            return vals.size == table.shape[0]

        def f_obs(vals, id_, md):
            return vals.size == table.shape[1]

        obs = table.filter(f_sample, inplace=False)
        self.assertEqual(obs, table)

        obs = table.filter(f_obs, 'observation', inplace=False)
        self.assertEqual(obs, table)

        def f(vals, id_, md):
            return (np.all(vals == [1, 0, 0]) or np.all(vals == [0, 0, 0]))

        obs = table.filter(f, inplace=False)
        exp = Table(np.array([[1, 0],
                              [0, 0],
                              [0, 0]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S4'])
        self.assertEqual(obs, exp)

        def f(vals, id_, md):
            return (np.all(vals == [0, 3, 4, 0]) or
                    np.all(vals == [0, 5, 0, 0]))

        obs = table.filter(f, 'observation', inplace=False)
        exp = Table(np.array([[0, 3, 4, 0],
                              [0, 5, 0, 0]]),
                    ['O1', 'O2'],
                    ['S1', 'S2', 'S3', 'S4'])
        self.assertNotEqual(obs, exp)

    def test_filter_id_state(self):
        def f(vals, id_, md):
            return id_[0] == 'b'

        filtered_table = self.st3.filter(f, inplace=False)
        filtered_table_2 = self.st3.filter(f, inplace=True)
        self.assertEqual(filtered_table._sample_index, {'b': 0})
        self.assertEqual(filtered_table._obs_index, {'2': 0, '3': 1})
        self.assertEqual(filtered_table_2._sample_index, {'b': 0})
        self.assertEqual(filtered_table_2._obs_index, {'2': 0, '3': 1})

    def test_filter_return_type(self):
        def f(vals, id_, md):
            return id_[0] == 'b'
        filtered_table = self.st3.filter(f, inplace=False)
        filtered_table_2 = self.st3.filter(f, inplace=True)
        self.assertEqual(filtered_table, filtered_table_2)
        self.assertTrue(filtered_table_2 is self.st3)

    def test_filter_general_sample(self):
        def f(vals, id_, md):
            return id_ == 'a'

        values = csr_matrix(np.array([[5.],
                                      [7.]]))
        exp_table = Table(values, ['1', '2'], ['a'],
                          [{'taxonomy': ['k__a', 'p__b']},
                           {'taxonomy': ['k__a', 'p__c']}],
                          [{'barcode': 'aatt'}])

        table = self.st_rich
        obs_table = table.filter(f, 'sample', inplace=False)
        self.assertEqual(obs_table, exp_table)

        def f_2(vals, id_, md):
            return np.all(vals == np.array([5, 7]))

        obs_table_2 = table.filter(f_2, 'sample', inplace=False)
        self.assertEqual(obs_table_2, exp_table)

    def test_filter_general_observation(self):
        def f(vals, id_, md):
            return md['taxonomy'][1] == 'p__c'

        values = csr_matrix(np.array([[7., 8.]]))
        exp_table = Table(values, ['2'], ['a', 'b'],
                          [{'taxonomy': ['k__a', 'p__c']}],
                          [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])
        table = self.st_rich
        obs_table = table.filter(f, 'observation', inplace=False)
        self.assertEqual(obs_table, exp_table)

        def f_2(vals, id_, md):
            return np.all(vals == np.array([7, 8]))
        obs_table_2 = table.filter(f_2, 'observation', inplace=False)
        self.assertEqual(obs_table_2, exp_table)

    def test_filter_sample_id(self):
        def f(vals, id_, md):
            return id_ == 'a'

        values = csr_matrix(np.array([[5.],
                                      [7.]]))
        exp_table = Table(values, ['1', '2'], ['a'],
                          [{'taxonomy': ['k__a', 'p__b']},
                           {'taxonomy': ['k__a', 'p__c']}],
                          [{'barcode': 'aatt'}])

        table = self.st_rich
        table.filter(f, 'sample')
        self.assertEqual(table, exp_table)

    def test_filter_sample_metadata(self):
        def f(vals, id_, md):
            return md['barcode'] == 'ttgg'
        values = csr_matrix(np.array([[6.],
                                      [8.]]))
        exp_table = Table(values, ['1', '2'], ['b'],
                          [{'taxonomy': ['k__a', 'p__b']},
                           {'taxonomy': ['k__a', 'p__c']}],
                          [{'barcode': 'ttgg'}])
        table = self.st_rich
        table.filter(f, 'sample')
        self.assertEqual(table, exp_table)

    def test_filter_sample_invert(self):
        def f(vals, id_, md):
            return md['barcode'] == 'aatt'
        values = csr_matrix(np.array([[6.],
                                      [8.]]))
        exp_table = Table(values, ['1', '2'], ['b'],
                          [{'taxonomy': ['k__a', 'p__b']},
                           {'taxonomy': ['k__a', 'p__c']}],
                          [{'barcode': 'ttgg'}])
        table = self.st_rich
        table.filter(f, 'sample', invert=True)
        self.assertEqual(table, exp_table)

    def test_filter_sample_remove_everything(self):
        def f(vals, id_, md):
            return False

        with errstate(empty='raise'), self.assertRaises(TableException):
            self.st_rich.filter(f, 'sample')

    def test_filter_observations_id(self):
        def f(vals, id_, md):
            return id_ == '1'

        values = csr_matrix(np.array([[5., 6.]]))
        exp_table = Table(values, ['1'], ['a', 'b'],
                          [{'taxonomy': ['k__a', 'p__b']}],
                          [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])
        table = self.st_rich
        table.filter(f, 'observation')
        self.assertEqual(table, exp_table)

    def test_filter_observations_metadata(self):
        def f(vals, id_, md):
            return md['taxonomy'][1] == 'p__c'

        values = csr_matrix(np.array([[7., 8.]]))
        exp_table = Table(values, ['2'], ['a', 'b'],
                          [{'taxonomy': ['k__a', 'p__c']}],
                          [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])
        table = self.st_rich
        table.filter(f, 'observation')
        self.assertEqual(table, exp_table)

    def test_filter_observations_invert(self):
        def f(vals, id_, md):
            return md['taxonomy'][1] == 'p__c'

        values = csr_matrix(np.array([[5., 6.]]))
        exp_table = Table(values, ['1'], ['a', 'b'],
                          [{'taxonomy': ['k__a', 'p__b']}],
                          [{'barcode': 'aatt'}, {'barcode': 'ttgg'}])
        table = self.st_rich
        table.filter(f, 'observation', invert=True)
        self.assertEqual(table, exp_table)

    def test_filter_observations_remove_everything(self):
        def f(vals, id_, md):
            return False

        with errstate(empty='raise'), self.assertRaises(TableException):
            self.st_rich.filter(f, 'observation')

    def test_subsample_edgecase_issue_952(self):
        # this file triggers an exception on Linux on subsample
        # with replacement where the pvals computed sum to > 1. It is a
        # subset of the data reported in issue 952, specifically constrained
        # to the first 10 features with any empty samples removed.
        path = 'test_data/edgecase_issue_952.biom'

        # ...existing logic for test_data, not ideal, but consistent
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        table = Table.from_hdf5(h5py.File(path, 'r'))
        os.chdir(cwd)

        obs = table.subsample(10, with_replacement=True)
        self.assertEqual(set(obs.sum('sample')), {10.0, })

    def test_subsample_same_seed_without_replacement(self):
        table = Table(np.array([[3, 1, 2], [0, 3, 4]]), ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        exp = table.subsample(2, seed=1234)
        for _ in range(100):
            obs = table.subsample(2, seed=1234)
            self.assertEqual(obs, exp)

    def test_subsample_same_seed_with_replacement(self):
        table = Table(np.array([[3, 1, 2], [0, 3, 4]]), ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        exp = table.subsample(2, seed=1234, with_replacement=True)
        for _ in range(100):
            obs = table.subsample(2, seed=1234, with_replacement=True)
            self.assertEqual(obs, exp)

    def test_subsample_by_id(self):
        table = Table(np.array([[3, 1, 2], [0, 3, 4]]), ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual_o1 = set()
        actual_o2 = set()
        for i in range(100):
            obs = table.subsample(2, by_id=True)
            actual_o1.add(tuple(obs.data('O1', 'observation')))
            actual_o2.add(tuple(obs.data('O2', 'observation')))
        self.assertEqual(actual_o1, {(3, 1), (1, 2), (3, 2)})
        self.assertEqual(actual_o2, {(0, 3), (3, 4), (0, 4)}),

    def test_subsample_by_id_observations_bug(self):
        table = Table(np.array([[3, 1, 2], [0, 3, 4]]).T, ['O1', 'O2', 'O3'],
                      ['S1', 'S2'])
        actual_o1 = set()
        actual_o2 = set()
        for i in range(100):
            obs = table.subsample(2, axis='observation', by_id=True)
            actual_o1.add(tuple(obs.data('S1', 'sample')))
            actual_o2.add(tuple(obs.data('S2', 'sample')))
        self.assertEqual(actual_o1, {(3, 1), (1, 2), (3, 2)})
        self.assertEqual(actual_o2, {(0, 3), (3, 4), (0, 4)}),

    def test_filter_using_list_of_ids(self):
        ids = ['S1', 'S4']
        obs = self.sparse_table.filter(ids, inplace=False)
        exp = Table(np.array([[1, 0],
                              [0, 0],
                              [0, 0]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S4'])
        self.assertEqual(obs, exp)

        ids = ['O1', 'O2']
        obs = self.sparse_table.filter(ids, 'observation', invert=True,
                                       inplace=False)
        exp = Table(np.array([[0, 5, 0, 0]]),
                    ['O3'],
                    ['S1', 'S2', 'S3', 'S4'])
        self.assertEqual(obs, exp)

    def test_filter_out_full_table(self):
        t = Table(np.asarray([[1, 2, 3],
                              [4, 5, 6]]),
                  ['a', 'b'], ['c', 'd', 'e'])
        t_sample = t.filter(ids_to_keep=[], axis='sample', inplace=False)
        t_obs = t.filter(ids_to_keep=[], axis='observation', inplace=False)

        self.assertEqual(t_sample.shape, (2, 0))
        self.assertEqual(t_obs.shape, (0, 3))

    def test_subsample(self):
        table = Table(np.array([[0, 5, 0]]), ['O1'], ['S1', 'S2', 'S3'])

        obs = table.subsample(5, axis='observation')
        npt.assert_equal(obs.data('O1', 'observation'), np.array([5]))
        self.assertEqual(obs.ids(), ['S2'])

        table = Table(np.array([[3, 1, 1], [0, 3, 3]]), ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual_o1 = set()
        actual_o2 = set()
        for i in range(100):
            obs = table.subsample(3)
            actual_o1.add(tuple(obs.data('O1', 'observation')))
            actual_o2.add(tuple(obs.data('O2', 'observation')))
        self.assertEqual(actual_o1, {(3, 0, 0), (3, 1, 0), (3, 0, 1),
                                     (3, 1, 1)})
        self.assertEqual(actual_o2, {(0, 3, 3), (0, 2, 3), (0, 3, 2),
                                     (0, 2, 2)})

    def test_subsample_md_copy_bug(self):
        """subsample would except when if metadata were present"""
        table = Table(np.array([[5, 5, 5]]), ['O1'], ['S1', 'S2', 'S3'],
                      [{'foo': 'bar'}], [{1: 2}, {3: 4}, {5: 6}])
        exp = table.copy()
        obs = table.subsample(5)
        self.assertEqual(obs, exp)

    def test_subsample_byid_with_replacement(self):
        dt = Table(np.array([[0, 1, 2], [2, 0, 1], [1, 2, 0]]),
                   ['O1', 'O2', 'O3'],
                   ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            dt.subsample(20, by_id=True, with_replacement=True)

    def test_subsample_without_replacement_unique_results(self):
        """
        As in scikit-bio. Given a vector of observations, the total number of
        unique subsamplings when n == sum(vector) - 1 should be equal to the
        number of unique categories of observations when the vector is
        subsampled without replacement. If the vector is subsamples *with*
        replacement, however, there should be more than 10 different
        possible subsamplings.
        """
        a = np.array([[2, 1, 2, 1, 8, 6, 3, 3, 5, 5], ]).T
        dt = Table(data=a, sample_ids=['S1', ],
                   observation_ids=[f'OTU{i:02d}' for i in range(10)])
        actual = set()
        for i in range(1000):
            obs = dt.subsample(35)
            actual.add(tuple(obs.data('S1')))
        self.assertEqual(len(actual), 10)

    def test_subsample_with_replacement_unique_results(self):
        """
        As in scikit-bio. Given a vector of observations, the total number of
        unique subsamplings when n == sum(vector) - 1 should be equal to the
        number of unique categories of observations when the vector is
        subsampled without replacement. If the vector is subsamples *with*
        replacement, however, there should be more than 10 different
        possible subsamplings.
        """
        a = np.array([[2, 1, 2, 1, 8, 6, 3, 3, 5, 5], ]).T
        dt = Table(data=a, sample_ids=['S1', ],
                   observation_ids=[f'OTU{i:02d}' for i in range(10)])
        actual = set()
        for i in range(1000):
            obs = dt.subsample(35, with_replacement=True)
            actual.add(tuple(obs.data('S1')))
        self.assertGreater(len(actual), 10)

    def test_subsample_with_replacement_n(self):
        dt = Table(np.array([[0, 1, 2], [2, 0, 1], [1, 2, 0]]),
                   ['O1', 'O2', 'O3'],
                   ['S1', 'S2', 'S3'])
        new_dt = dt.subsample(20, with_replacement=True)
        new_counts = np.unique(new_dt.sum('sample'))
        self.assertEqual(new_counts.shape[0], 1)
        self.assertEqual(new_counts[0], 20)

    def test_pa(self):
        exp = Table(np.array([[1, 1], [1, 0]]), ['5', '6'], ['a', 'b'])
        self.st7.pa()
        self.assertEqual(self.st7, exp)

    def test_pa_with_neg(self):
        t = Table(np.array([[-10, 7], [0, -0.1]]), ['5', '6'], ['a', 'b'])
        exp = Table(np.array([[1, 1], [0, 1]]), ['5', '6'], ['a', 'b'])
        t.pa()
        self.assertEqual(t, exp)

    def test_pa_works_if_something_has_been_zeroed(self):
        exp = Table(np.array([[0, 1], [1, 0]]), ['5', '6'], ['a', 'b'])
        self.st7._data[0, 0] = 0
        self.st7.pa()
        self.assertEqual(self.st7, exp)

    def test_transform_return_type(self):
        def f(data, id_, md):
            return data / 2.

        filtered_table = self.st3.transform(f, inplace=False)
        filtered_table_2 = self.st3.transform(f, inplace=True)
        self.assertEqual(filtered_table, filtered_table_2)
        self.assertTrue(filtered_table_2 is self.st3)

    def test_transform_observation(self):
        """Transform axis by arbitrary function"""
        # Transform observations by arbitrary function
        def obs_transform_f(v, id, md):
            return np.where(v >= 7, 1, 0)
        sp_sd = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 1}
        exp = Table(sp_sd, ['1', '2'], ['a', 'b'])
        self.st1.transform(obs_transform_f, axis='observation')
        self.assertEqual(self.st1, exp)

    def test_transform_sample(self):
        # """Transform samples by arbitrary function"""
        def sample_transform_f(v, id, md):
            return np.where(v >= 6, 1, 0)

        sp_sd = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 1}
        exp = Table(sp_sd, ['1', '2'], ['a', 'b'])
        self.st1.transform(sample_transform_f)
        self.assertEqual(self.st1, exp)

        # Raises UnknownAxisError if a invalid axis is passed
        with self.assertRaises(UnknownAxisError):
            self.st1.transform(sample_transform_f, axis='foo')

    def test_rank_observation_by_sample(self):
        """rank observations by sample"""
        data = np.array([[99, 12, 8],
                         [0, 42, 7],
                         [112, 42, 6],
                         [5, 75, 5]])
        data_exp = np.array([[2., 1., 4.],
                             [0., 2.5, 3.],
                             [3., 2.5, 2.],
                             [1., 4., 1.]])

        st = Table(data, sample_ids=['s1', 's2', 's3'],
                   observation_ids=['o1', 'o2', 'o3', 'o4'])
        exp = Table(data_exp, sample_ids=['s1', 's2', 's3'],
                    observation_ids=['o1', 'o2', 'o3', 'o4'])
        st.rankdata(axis='sample')
        self.assertEqual(st, exp)

    def test_rank_observation_by_sample_alt_method(self):
        """rank observations by sample with alt method"""
        data = np.array([[99, 12, 8],
                         [0, 42, 7],
                         [112, 42, 6],
                         [5, 75, 5]])
        data_exp = np.array([[2., 1., 4.],
                             [0., 2., 3.],
                             [3., 2., 2.],
                             [1., 4., 1.]])

        st = Table(data, sample_ids=['s1', 's2', 's3'],
                   observation_ids=['o1', 'o2', 'o3', 'o4'])
        exp = Table(data_exp, sample_ids=['s1', 's2', 's3'],
                    observation_ids=['o1', 'o2', 'o3', 'o4'])
        st.rankdata(axis='sample', method='min')
        self.assertEqual(st, exp)

    def test_rank_sample_by_observation(self):
        """rank samples by observation"""
        data = np.array([[99, 12, 8],
                         [0, 42, 7],
                         [112, 42, 6],
                         [5, 75, 5]])
        data_exp = np.array([[3., 2., 1.],
                             [0., 2., 1.],
                             [3., 2., 1.],
                             [1.5, 3., 1.5]])

        st = Table(data, sample_ids=['s1', 's2', 's3'],
                   observation_ids=['o1', 'o2', 'o3', 'o4'])
        exp = Table(data_exp, sample_ids=['s1', 's2', 's3'],
                    observation_ids=['o1', 'o2', 'o3', 'o4'])
        st.rankdata(axis='observation')
        self.assertEqual(st, exp)

    def test_norm_observation_by_sample(self):
        """normalize observations by sample"""
        data = {(0, 0): 2, (0, 1): 0, (1, 0): 6, (1, 1): 1}
        data_exp = {(0, 0): 0.25, (0, 1): 0.0, (1, 0): 0.75, (1, 1): 1.0}

        st = Table(data, ['1', '2'], ['a', 'b'])
        exp = Table(data_exp, ['1', '2'], ['a', 'b'])
        st.norm()
        self.assertEqual(st, exp)

    def test_norm_sample_by_observation(self):
        """normalize sample by observation"""
        data = {(0, 0): 0, (0, 1): 2, (1, 0): 2, (1, 1): 6}
        data_exp = {(0, 0): 0.0, (0, 1): 1.0, (1, 0): 0.25, (1, 1): 0.75}
        st = Table(data, ['1', '2'], ['a', 'b'])
        exp = Table(data_exp, ['1', '2'], ['a', 'b'])
        st.norm(axis='observation')
        self.assertEqual(st, exp)

    def test_collapse_observations_by_metadata_one_to_many_strict(self):
        """Collapse observations by arbitary metadata"""
        dt_rich = Table(np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
                        ['1', '2', '3'], ['a', 'b', 'c'],
                        [{'pathways': [['a', 'bx'], ['a', 'd']]},
                         {'pathways': [['a', 'bx'], ['a', 'c']]},
                         {'pathways': [['a']]}],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}])
        exp_cat2 = Table(np.array([[13, 15, 17], [8, 9, 10], [5, 6, 7]]),
                         ['bx', 'c', 'd'], ['a', 'b', 'c'],
                         [{'Path': ['a', 'bx']},
                          {'Path': ['a', 'c']},
                          {'Path': ['a', 'd']}],
                         [{'barcode': 'aatt'},
                          {'barcode': 'ttgg'},
                          {'barcode': 'aatt'}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo, foo[1])

        obs_cat2 = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1, one_to_many=True,
            strict=False, axis='observation').sort(axis='observation')
        self.assertEqual(obs_cat2, exp_cat2)

        with self.assertRaises(IndexError):
            dt_rich.collapse(
                bin_f, norm=False, min_group_size=1, one_to_many=True,
                strict=True, axis='observation')

    def test_collapse_observations_by_metadata_one_to_many(self):
        """Collapse observations by arbitary metadata"""
        dt_rich = Table(np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13],
                                  [14, 15, 16]]),
                        ['1', '2', '3', '4'], ['a', 'b', 'c'],
                        [{'pathways': [['a', 'bx'], ['a', 'd']]},
                         {'pathways': [['a', 'bx'], ['a', 'c']]},
                         {'pathways': [['a', 'c']]},
                         {'pathways': [['a', 'c']]}],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}])
        exp_cat2 = Table(np.array([[13, 15, 17], [33, 36, 39], [5, 6, 7]]),
                         ['bx', 'c', 'd'], ['a', 'b', 'c'],
                         [{'Path': ['a', 'bx']},
                          {'Path': ['a', 'c']},
                          {'Path': ['a', 'd']}],
                         [{'barcode': 'aatt'},
                          {'barcode': 'ttgg'},
                          {'barcode': 'aatt'}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo, foo[-1])

        obs_cat2 = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            one_to_many=True, axis='observation').sort(axis='observation')
        self.assertEqual(obs_cat2, exp_cat2)

        dt_rich = Table(np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
                        ['1', '2', '3'], ['a', 'b', 'c'],
                        [{'pathways': [['a', 'b'], ['a', 'd']]},
                         {'pathways': [['a', 'b'], ['a', 'c']]},
                         {'pathways': [['a', 'c']]}],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}])
        exp_cat1 = Table(np.array([[37, 42, 47]]),
                         ['a'], ['a', 'b', 'c'],
                         [{'Path': ['a']}],
                         [{'barcode': 'aatt'},
                          {'barcode': 'ttgg'},
                          {'barcode': 'aatt'}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo[:1], foo[0])

        obs_cat1 = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            one_to_many=True, axis='observation').sort(axis='observation')
        self.assertEqual(obs_cat1, exp_cat1)

        # Test out include_collapsed_metadata=False.
        exp = Table(np.array([[37, 42, 47]]),
                    ['a'], ['a', 'b', 'c'], None,
                    [{'barcode': 'aatt'},
                     {'barcode': 'ttgg'},
                     {'barcode': 'aatt'}])
        obs = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1, one_to_many=True,
            include_collapsed_metadata=False,
            axis='observation').sort(axis='observation')
        self.assertEqual(obs, exp)

        # Test out constructor.
        obs = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1, one_to_many=True,
            include_collapsed_metadata=False,
            axis='observation').sort(axis='observation')
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), Table)

    def test_collapse_observations_by_metadata_one_to_many_divide(self):
        """Collapse observations by 1-M metadata using divide mode"""
        dt_rich = Table(np.array([[1, 6, 7], [8, 0, 10], [11, 12, 13]]),
                        ['1', '2', '3'],
                        ['a', 'b', 'c'],
                        [{'pathways': [['a', 'bx'], ['a', 'd']]},
                         {'pathways': [['a', 'bx'], ['a', 'c']]},
                         {'pathways': [['a', 'c']]}],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}])
        exp = Table(np.array([[4.5, 3, 8.5], [15, 12, 18], [0.5, 3, 3.5]]),
                    ['bx', 'c', 'd'],
                    ['a', 'b', 'c'],
                    [{'Path': ['a', 'bx']},
                     {'Path': ['a', 'c']},
                     {'Path': ['a', 'd']}],
                    [{'barcode': 'aatt'},
                     {'barcode': 'ttgg'},
                     {'barcode': 'aatt'}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo, foo[-1])

        obs = dt_rich.collapse(
            bin_f, norm=False, one_to_many=True,
            one_to_many_mode='divide',
            axis='observation').sort(axis='observation')
        self.assertEqual(obs, exp)

        # Test skipping some observation metadata (strict=False).
        dt_rich = Table(
            np.array([[5.0, 6.0, 7], [8, 9, 10], [11, 12, 13.0]]),
            ['1', '2', '3'], ['a', 'b', 'c'],
            [{'pathways': [['a', 'bx'], ['a', 'd']]},
             {'pathways': [['a', 'bx'], ['a', 'c'], ['z']]},
             {'pathways': [['a']]}],
            [{'barcode': 'aatt'},
             {'barcode': 'ttgg'},
             {'barcode': 'aatt'}])
        exp = Table(np.array([[6.5, 7.5, 8.5], [4, 4.5, 5], [2.5, 3, 3.5]]),
                    ['bx', 'c', 'd'], ['a', 'b', 'c'],
                    [{'Path': ['a', 'bx']},
                     {'Path': ['a', 'c']},
                     {'Path': ['a', 'd']}],
                    [{'barcode': 'aatt'},
                     {'barcode': 'ttgg'},
                     {'barcode': 'aatt'}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo, foo[1])

        obs = dt_rich.collapse(
            bin_f, norm=False, one_to_many=True, one_to_many_mode='divide',
            strict=False, axis='observation').sort(axis='observation')

        self.assertEqual(obs, exp)

        with self.assertRaises(IndexError):
            dt_rich.collapse(
                bin_f, norm=False, one_to_many=True, one_to_many_mode='divide',
                strict=True, axis='observation')

        # Invalid one_to_many_mode.
        with self.assertRaises(ValueError):
            dt_rich.collapse(
                bin_f, norm=False, one_to_many=True, one_to_many_mode='foo',
                axis='observation')

    def test_collapse_median(self):
        table = Table(
            np.array([[5, 6, 7],
                      [1, 2, 3],
                      [8, 9, 10],
                      [1, 2.5, 1],
                      [11, 12, 13],
                      [2, 3, 10]]),
            ['a', 'b', 'c', 'd', 'e', 'f'],
            ['s1', 's2', 's3'])

        # two partitions, (a, c, e) and (b, d, f)
        def partition_f(id_, md):
            return id_ in {'b', 'd', 'f'}

        def collapse_f(t, axis):
            return np.asarray([np.median(v) for v in t.iter_data(dense=True)])

        obs = table.collapse(partition_f, collapse_f, axis='observation',
                             norm=False)
        exp = Table(np.array([[8, 9, 10], [1, 2.5, 3]]),
                    [False, True],
                    ['s1', 's2', 's3'],
                    [{'collapsed_ids': ['a', 'c', 'e']},
                     {'collapsed_ids': ['b', 'd', 'f']}])
        self.assertEqual(obs, exp)

    def test_collapse_observations_by_metadata(self):
        """Collapse observations by arbitrary metadata"""
        dt_rich = Table(
            np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
            ['1', '2', '3'], ['a', 'b', 'c'],
            [{'taxonomy': ['k__a', 'p__b']},
             {'taxonomy': ['k__a', 'p__c']},
             {'taxonomy': ['k__a', 'p__c']}],
            [{'barcode': 'aatt'},
             {'barcode': 'ttgg'},
             {'barcode': 'aatt'}])
        exp_phy = Table(np.array([[5, 6, 7], [19, 21, 23]]),
                        ['p__b', 'p__c'], ['a', 'b', 'c'],
                        [{'collapsed_ids': ['1']},
                         {'collapsed_ids': ['2', '3']}],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}])

        def bin_f(id_, x):
            return x['taxonomy'][1]

        obs_phy = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            axis='observation').sort(axis='observation')
        self.assertEqual(obs_phy, exp_phy)

        exp_king = Table(np.array([[24, 27, 30]]),
                         ['k__a'], ['a', 'b', 'c'],
                         [{'collapsed_ids': ['1', '2', '3']}],
                         [{'barcode': 'aatt'},
                          {'barcode': 'ttgg'},
                          {'barcode': 'aatt'}])

        def bin_f(id_, x):
            return x['taxonomy'][0]

        obs_king = dt_rich.collapse(bin_f, norm=False, axis='observation')
        self.assertEqual(obs_king, exp_king)

        with errstate(all='raise'), self.assertRaises(TableException):
            dt_rich.collapse(bin_f, min_group_size=10, axis='observation')

        # Test out include_collapsed_metadata=False.
        exp = Table(np.array([[24, 27, 30]]),
                    ['k__a'],
                    ['a', 'b', 'c'], None,
                    [{'barcode': 'aatt'},
                     {'barcode': 'ttgg'},
                     {'barcode': 'aatt'}])
        obs = dt_rich.collapse(bin_f, norm=False,
                               include_collapsed_metadata=False,
                               axis='observation')
        self.assertEqual(obs, exp)

        # Test out constructor.
        obs = dt_rich.collapse(bin_f, norm=False,
                               include_collapsed_metadata=False,
                               axis='observation')
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), Table)

    def test_collapse_samples_by_metadata(self):
        """Collapse samples by arbitrary metadata"""
        dt_rich = Table(
            np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
            ['1', '2', '3'], ['a', 'b', 'c'],
            [{'taxonomy': ['k__a', 'p__b']},
             {'taxonomy': ['k__a', 'p__c']},
             {'taxonomy': ['k__a', 'p__c']}],
            [{'barcode': 'aatt'},
             {'barcode': 'ttgg'},
             {'barcode': 'aatt'}])
        exp_bc = Table(
            np.array([[12, 6], [18, 9], [24, 12]]),
            ['1', '2', '3'], ['aatt', 'ttgg'],
            [{'taxonomy': ['k__a', 'p__b']},
             {'taxonomy': ['k__a', 'p__c']},
             {'taxonomy': ['k__a', 'p__c']}],
            [{'collapsed_ids': ['a', 'c']},
             {'collapsed_ids': ['b']}])

        def bin_f(id_, x):
            return x['barcode']

        obs_bc = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            axis='sample').sort(axis='sample')
        self.assertEqual(obs_bc, exp_bc)

        with errstate(all='raise'), self.assertRaises(TableException):
            dt_rich.collapse(bin_f, min_group_size=10)

        # Test out include_collapsed_metadata=False.
        exp = Table(np.array([[12, 6], [18, 9], [24, 12]]),
                    ['1', '2', '3'],
                    ['aatt', 'ttgg'],
                    [{'taxonomy': ['k__a', 'p__b']},
                     {'taxonomy': ['k__a', 'p__c']},
                     {'taxonomy': ['k__a', 'p__c']}],
                    None)

        obs = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            include_collapsed_metadata=False).sort(axis='sample')
        self.assertEqual(obs, exp)

        # Test out constructor.
        obs = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            include_collapsed_metadata=False).sort(axis='sample')
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), Table)

    def test_collapse_samples_by_metadata_one_to_many_strict(self):
        """Collapse samples by arbitary metadata"""
        dt_rich = Table(np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
                        ['1', '2', '3'],
                        ['XXa', 'XXb', 'XXc'],
                        [{'other': 'aatt'},
                         {'other': 'ttgg'},
                         {'other': 'aatt'}],
                        [{'foo': [['a', 'b'], ['a', 'd']]},
                         {'foo': [['a', 'b'], ['a', 'c']]},
                         {'foo': [['a']]}])
        exp_cat2 = Table(np.array([[11, 17, 23], [6, 9, 12], [5, 8, 11]]).T,
                         ['1', '2', '3'],
                         ['b', 'c', 'd'],
                         [{'other': 'aatt'},
                          {'other': 'ttgg'},
                          {'other': 'aatt'}],
                         [{'Path': ['a', 'b']},
                          {'Path': ['a', 'c']},
                          {'Path': ['a', 'd']}])

        def bin_f(id_, x):
            for foo in x['foo']:
                yield (foo, foo[1])

        obs_cat2 = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1, one_to_many=True,
            strict=False).sort(axis='observation')
        self.assertEqual(obs_cat2, exp_cat2)

        self.assertRaises(IndexError, dt_rich.collapse, bin_f,
                          norm=False, min_group_size=1, one_to_many=True,
                          strict=True)

    def test_collapse_samples_by_metadata_one_to_many_divide(self):
        """Collapse samples by 1-M metadata using divide mode"""
        dt_rich = Table(np.array([[1, 8, 11], [6, 0, 12], [7, 10, 13]]),
                        ['a', 'b', 'c'],
                        ['1', '2', '3'],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}],
                        [{'pathways': [['a', 'bx'], ['a', 'd']]},
                         {'pathways': [['a', 'bx'], ['a', 'c']]},
                         {'pathways': [['a', 'c']]}])
        exp = Table(np.array([[4.5, 15, 0.5], [3, 12, 3], [8.5, 18, 3.5]]),
                    ['a', 'b', 'c'],
                    ['bx', 'c', 'd'],
                    [{'barcode': 'aatt'},
                     {'barcode': 'ttgg'},
                     {'barcode': 'aatt'}],
                    [{'Path': ['a', 'bx']},
                     {'Path': ['a', 'c']},
                     {'Path': ['a', 'd']}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo, foo[-1])

        obs = dt_rich.collapse(
            bin_f, norm=False, one_to_many=True,
            one_to_many_mode='divide').sort(axis='sample')
        self.assertEqual(obs, exp)

        # Test skipping some sample metadata (strict=False).
        dt_rich = Table(np.array([[5.0, 8, 11], [6.0, 9, 12], [7, 10, 13.0]]),
                        ['a', 'b', 'c'],
                        ['1', '2', '3'],
                        [{'barcode': 'aatt'},
                         {'barcode': 'ttgg'},
                         {'barcode': 'aatt'}],
                        [{'pathways': [['a', 'bx'], ['a', 'd']]},
                         {'pathways': [['a', 'bx'], ['a', 'c'], ['z']]},
                         {'pathways': [['a']]}])
        exp = Table(np.array([[6.5, 4, 2.5], [7.5, 4.5, 3], [8.5, 5, 3.5]]),
                    ['a', 'b', 'c'],
                    ['bx', 'c', 'd'],
                    [{'barcode': 'aatt'},
                     {'barcode': 'ttgg'},
                     {'barcode': 'aatt'}],
                    [{'Path': ['a', 'bx']},
                     {'Path': ['a', 'c']},
                     {'Path': ['a', 'd']}])

        def bin_f(id_, x):
            for foo in x['pathways']:
                yield (foo, foo[1])

        obs = dt_rich.collapse(
            bin_f, norm=False, one_to_many=True, one_to_many_mode='divide',
            strict=False).sort(axis='sample')

        self.assertEqual(obs, exp)

        with self.assertRaises(IndexError):
            dt_rich.collapse(bin_f, norm=False,
                             one_to_many=True,
                             one_to_many_mode='divide',
                             strict=True)

        # Invalid one_to_many_mode.
        with self.assertRaises(ValueError):
            dt_rich.collapse(bin_f, norm=False,
                             one_to_many=True,
                             one_to_many_mode='foo')

    def test_collapse_samples_by_metadata_one_to_many(self):
        """Collapse samples by arbitary metadata"""
        dt_rich = Table(np.array([[5, 6, 7],
                                  [8, 9, 10],
                                  [11, 12, 13]]),
                        ['1', '2', '3'],
                        ['XXa', 'XXb', 'XXc'],
                        [{'other': 'aatt'},
                         {'other': 'ttgg'},
                         {'other': 'aatt'}],
                        [{'foo': [['a', 'b'], ['a', 'd']]},
                         {'foo': [['a', 'b'], ['a', 'c']]},
                         {'foo': [['a', 'c']]}])
        exp_cat2 = Table(
            np.array([[11, 17, 23], [13, 19, 25], [5, 8, 11]]).T,
            ['1', '2', '3'],
            ['b', 'c', 'd'],
            [{'other': 'aatt'},
             {'other': 'ttgg'},
             {'other': 'aatt'}],
            [{'Path': ['a', 'b']},
             {'Path': ['a', 'c']},
             {'Path': ['a', 'd']}])

        def bin_f(id_, x):
            for foo in x['foo']:
                yield (foo, foo[-1])

        obs_cat2 = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            one_to_many=True, axis='sample').sort(axis='observation')

        self.assertEqual(obs_cat2, exp_cat2)

        dt_rich = Table(
            np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
            ['1', '2', '3'], ['a', 'b', 'c'],
            [{'other': 'aatt'},
             {'other': 'ttgg'},
             {'other': 'aatt'}],
            [{'foo': [['a', 'b'], ['a', 'd']]},
             {'foo': [['a', 'b'], ['a', 'c']]},
             {'foo': [['a', 'c']]}])
        exp_cat1 = Table(np.array([[29, 44, 59]]).T,
                         ['1', '2', '3'], ['a'],
                         [{'other': 'aatt'},
                          {'other': 'ttgg'},
                          {'other': 'aatt'}],
                         [{'Path': ['a']}])

        def bin_f(id_, x):
            for foo in x['foo']:
                yield (foo[:1], foo[0])

        obs_cat1 = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1,
            one_to_many=True, axis='sample').sort(axis='observation')
        self.assertEqual(obs_cat1, exp_cat1)

        # Test out include_collapsed_metadata=False.
        exp = Table(np.array([[29, 44, 59]]).T,
                    ['1', '2', '3'],
                    ['a'],
                    [{'other': 'aatt'},
                     {'other': 'ttgg'},
                     {'other': 'aatt'}],
                    None)
        obs = dt_rich.collapse(
            bin_f, norm=False, min_group_size=1, one_to_many=True,
            include_collapsed_metadata=False,
            axis='sample').sort(axis='observation')
        self.assertEqual(obs, exp)

        # Test out constructor.
        obs = dt_rich.collapse(bin_f, norm=False, min_group_size=1,
                               one_to_many=True,
                               include_collapsed_metadata=False,
                               axis='sample').sort(axis='observation')
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), Table)

    def test_from_json_issue_697(self):
        t = Table({}, [], [])
        serialized = t.to_json('foo')
        reloaded = Table.from_json(loads(serialized))
        self.assertEqual(t, reloaded)
        self.assertEqual(reloaded.generated_by, 'foo')

    def test_to_json_int64_metadata_issue_886(self):
        t = example_table.copy()
        t.add_metadata({'S1': {'X': np.int64(1)},
                        'S2': {'X': np.int64(2)},
                        'S3': {'X': np.int64(3)}})
        t.to_json('test')

    def test_to_json_empty(self):
        t = Table({}, [], [])
        serialized = t.to_json('foo')
        reloaded = Table.from_json(loads(serialized))
        self.assertEqual(t, reloaded)

    def test_to_json_dense_int(self):
        """Get a BIOM format string for a dense table of integers"""
        # check by round trip
        obs_ids = list(map(str, range(5)))
        samp_ids = list(map(str, range(10)))
        obs_md = [{'foo': i} for i in range(5)]
        samp_md = [{'bar': i} for i in range(10)]
        data = np.reshape(np.arange(50), (5, 10))

        # using Table type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md)

        # verify that we can parse still
        t2 = parse_biom_table(StringIO(t.to_json('asd')))

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_dense_float(self):
        """Get a BIOM format string for a dense table of floats"""
        # check by round trip
        obs_ids = ['a', 'b']
        samp_ids = ['c', 'd']
        obs_md = [{'foo': i} for i in range(2)]
        samp_md = [{'bar': i} for i in range(2)]
        data = np.array([[0.01, 1.5], [0.0, 0.79]])

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md)

        # verify that we can parse still
        t2 = parse_biom_table(StringIO(t.to_json('asd')))

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_dense_int_directio(self):
        """Get a BIOM format string for a dense table of integers"""
        # check by round trip
        obs_ids = list(map(str, range(5)))
        samp_ids = list(map(str, range(10)))
        obs_md = [{'foo': i} for i in range(5)]
        samp_md = [{'bar': i} for i in range(10)]
        data = np.reshape(np.arange(50), (5, 10))

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md)

        # verify that we can parse still
        io = StringIO()
        t.to_json('asd', direct_io=io)
        io.seek(0)
        t2 = parse_biom_table(io)

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_dense_float_directio(self):
        """Get a BIOM format string for a dense table of floats"""
        # check by round trip
        obs_ids = ['a', 'b']
        samp_ids = ['c', 'd']
        obs_md = [{'foo': i} for i in range(2)]
        samp_md = [{'bar': i} for i in range(2)]
        data = np.array([[0.01, 1.5], [0.0, 0.79]])

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md)

        # verify that we can parse still
        io = StringIO()
        t.to_json('asd', direct_io=io)
        io.seek(0)
        t2 = parse_biom_table(io)

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_sparse_int(self):
        """Get a BIOM format string for a sparse table of integers"""
        # check by round trip
        obs_ids = list(map(str, range(5)))
        samp_ids = list(map(str, range(10)))
        obs_md = [{'foo': i} for i in range(5)]
        samp_md = [{'bar': i} for i in range(10)]
        data = [[0, 0, 10], [1, 1, 11], [2, 2, 12], [3, 3, 13], [4, 4, 14],
                [3, 5, 15], [2, 6, 16], [1, 7, 18], [0, 8, 19], [1, 9, 20]]

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md, obs_md)

        # verify that we can parse still
        t2 = parse_biom_table(StringIO(t.to_json('asd')))

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_sparse_float(self):
        """Get a BIOM format string for a sparse table of floats"""
        # check by round trip
        obs_ids = ['a', 'b']
        samp_ids = ['c', 'd']
        obs_md = [{'foo': i} for i in range(2)]
        samp_md = [{'bar': i} for i in range(2)]
        data = [[0, 0, 0.01], [0, 1, 1.5], [1, 0, 0.0], [1, 1, 0.79]]

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md, obs_md)

        # verify that we can parse still
        t2 = parse_biom_table(StringIO(t.to_json('asd')))

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_sparse_float_creation_date(self):
        """Verify we can inject a creation date"""
        # check by round trip
        obs_ids = ['a', 'b']
        samp_ids = ['c', 'd']
        obs_md = [{'foo': i} for i in range(2)]
        samp_md = [{'bar': i} for i in range(2)]
        data = [[0, 0, 0.01], [0, 1, 1.5], [1, 0, 0.0], [1, 1, 0.79]]
        current = datetime.now()

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md, obs_md)

        # verify that we can parse still
        t2 = parse_biom_table(StringIO(t.to_json('asd',
                                                 creation_date=current)))

        # verify that the tables are the same
        self.assertEqual(t, t2)
        self.assertEqual(t2.create_date, current)

    def test_to_json_sparse_int_directio(self):
        """Get a BIOM format string for a sparse table of integers"""
        # check by round trip
        obs_ids = list(map(str, range(5)))
        samp_ids = list(map(str, range(10)))
        obs_md = [{'foo': i} for i in range(5)]
        samp_md = [{'bar': i} for i in range(10)]
        data = [[0, 0, 10], [1, 1, 11], [2, 2, 12], [3, 3, 13], [4, 4, 14],
                [3, 5, 15], [2, 6, 16], [1, 7, 18], [0, 8, 19], [1, 9, 20]]

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md, obs_md)

        # verify that we can parse still
        io = StringIO()
        t.to_json('asd', direct_io=io)
        io.seek(0)
        t2 = parse_biom_table(io)

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_to_json_sparse_float_directio(self):
        """Get a BIOM format string for a sparse table of floats"""
        # check by round trip
        obs_ids = ['a', 'b']
        samp_ids = ['c', 'd']
        obs_md = [{'foo': i} for i in range(2)]
        samp_md = [{'bar': i} for i in range(2)]
        data = [[0, 0, 0.01], [0, 1, 1.5], [1, 0, 0.0], [1, 1, 0.79]]

        # using OTUTable type to support parsing round trip
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md)

        # verify that we can parse still
        io = StringIO()
        t.to_json('asd', direct_io=io)
        io.seek(0)
        t2 = parse_biom_table(io)

        # verify that the tables are the same
        self.assertEqual(t, t2)

    def test_extract_data_from_tsv(self):
        """Parses a classic table

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        input = legacy_otu_table1.splitlines()
        samp_ids = ['Fing', 'Key', 'NA']
        obs_ids = ['0', '1', '7', '3', '4']
        metadata = [
            'Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; '
            'Propionibacterium',
            'Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillal'
            'es; Lactobacillales; Streptococcaceae; Streptococcus',
            'Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Coryneb'
            'acteriaceae',
            'Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococca'
            'ceae',
            'Bacteria; Cyanobacteria; Chloroplasts; vectors']
        md_name = 'Consensus Lineage'
        data = [[0, 0, 19111], [0, 1, 44536], [0, 2, 42],
                [1, 0, 1216], [1, 1, 3500], [1, 2, 6],
                [2, 0, 1803], [2, 1, 1184], [2, 2, 2],
                [3, 0, 1722], [3, 1, 4903], [3, 2, 17],
                [4, 0, 589], [4, 1, 2074], [4, 2, 34]]

        exp = (samp_ids, obs_ids, data, metadata, md_name)
        obs = Table._extract_data_from_tsv(input, dtype=int)
        npt.assert_equal(obs, exp)

    def test_extract_data_from_tsv_bad_metadata(self):
        input = legacy_otu_table_bad_metadata.splitlines()
        samp_ids = ['Fing', 'Key', 'NA']
        obs_ids = ['0', '1', '7', '3', '4']
        metadata = [
            '',
            'Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillal'
            'es; Lactobacillales; Streptococcaceae; Streptococcus',
            'Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Coryneb'
            'acteriaceae',
            'Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococca'
            'ceae',
            'Bacteria; Cyanobacteria; Chloroplasts; vectors']
        md_name = 'Consensus Lineage'
        data = [[0, 0, 19111], [0, 1, 44536], [0, 2, 42],
                [1, 0, 1216], [1, 1, 3500], [1, 2, 6],
                [2, 0, 1803], [2, 1, 1184], [2, 2, 2],
                [3, 0, 1722], [3, 1, 4903], [3, 2, 17],
                [4, 0, 589], [4, 1, 2074], [4, 2, 34]]

        exp = (samp_ids, obs_ids, data, metadata, md_name)
        obs = Table._extract_data_from_tsv(input, dtype=int)
        npt.assert_equal(obs, exp)

        # and assert the exact identified bug in #827 is resolved
        input = extract_tsv_bug.splitlines()
        samp_ids = ['s1', 's2']
        obs_ids = ['1', '2', '3']
        metadata = [
            '',
            'k__test;p__test',
            'k__test;p__test']
        md_name = 'taxonomy'
        data = [[0, 0, 123], [0, 1, 32],
                [1, 0, 315], [1, 1, 3],
                [2, 1, 22]]

        exp = (samp_ids, obs_ids, data, metadata, md_name)
        obs = Table._extract_data_from_tsv(input, dtype=int)
        npt.assert_equal(obs, exp)

    def test_identify_bad_value(self):
        pos = [str(i) for i in range(10)]
        exp = (None, None)
        obs = _identify_bad_value(int, pos)
        self.assertEqual(obs, exp)

        neg = list('01234x6789')
        exp = ('x', 5)
        obs = _identify_bad_value(int, neg)
        self.assertEqual(obs, exp)

    def test_extract_data_from_tsv_badvalue_complaint(self):
        tsv = ['#OTU ID\ta\tb', '1\t2\t3', '2\tfoo\t6']

        msg = "Invalid value on line 2, column 1, value foo"
        with self.assertRaisesRegex(TypeError, msg):
            Table._extract_data_from_tsv(tsv, dtype=int)

    def test_partition_remove_empty(self):
        t = Table(np.array([[0, 1, 2],
                            [3, 0, 0],
                            [4, 0, 0]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        part_f = lambda i, m: i == 'S1'  # noqa
        obs = dict(t.partition(part_f, remove_empty=True))
        exp = {True: Table(np.array([[3, ], [4, ]]), ['O2', 'O3'], ['S1', ]),
               False: Table(np.array([[1, 2]]), ['O1', ], ['S2', 'S3'])}
        self.assertEqual(obs, exp)

    def test_partition_ignore_none_true(self):
        t = Table(np.array([[0, 1, 2],
                            [3, 0, 0],
                            [4, 0, 0]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        part_f = lambda i, m: True if i == 'S1' else None  # noqa
        obs = dict(t.partition(part_f, ignore_none=True))
        exp = {True: Table(np.array([[0, ], [3, ], [4, ]]),
                           ['O1', 'O2', 'O3'], ['S1', ])}
        self.assertEqual(obs, exp)

    def test_partition_ignore_none_false(self):
        t = Table(np.array([[0, 1, 2],
                            [3, 0, 0],
                            [4, 0, 0]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        part_f = lambda i, m: True if i == 'S1' else None  # noqa
        obs = dict(t.partition(part_f, ignore_none=False))
        exp = {True: Table(np.array([[0, ], [3, ], [4, ]]),
                           ['O1', 'O2', 'O3'], ['S1', ]),
               None: Table(np.array([[1, 2], [0, 0], [0, 0]]),
                           ['O1', 'O2', 'O3'], ['S2', 'S3'])}
        self.assertEqual(obs, exp)

    def test_partition_dict_ids_to_groups(self):
        t = Table(np.array([[0, 1, 2],
                            [3, 0, 0],
                            [4, 0, 0]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        by_dict = {'S1': 'foo',
                   'S2': 'bar',
                   'S3': 'foo'}
        exp = {'foo': Table(np.array([[0, 2], [3, 0], [4, 0]]),
                            ['O1', 'O2', 'O3'],
                            ['S1', 'S3']),
               'bar': Table(np.array([[1, ], [0, ], [0, ]]),
                            ['O1', 'O2', 'O3'],
                            ['S2', ])}
        obs = dict(t.partition(by_dict))
        self.assertEqual(obs, exp)

    def test_partition_dict_groups_to_ids(self):
        t = Table(np.array([[0, 1, 2],
                            [3, 0, 0],
                            [4, 0, 0]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        by_dict_group = {'foo': ['S1', 'S3'],
                         'bar': ['S2', ]}
        exp = {'foo': Table(np.array([[0, 2], [3, 0], [4, 0]]),
                            ['O1', 'O2', 'O3'],
                            ['S1', 'S3']),
               'bar': Table(np.array([[1, ], [0, ], [0, ]]),
                            ['O1', 'O2', 'O3'],
                            ['S2', ])}
        obs = dict(t.partition(by_dict_group))
        self.assertEqual(obs, exp)

    def test_bin_samples_by_metadata(self):
        """Yield tables binned by sample metadata"""
        def f(id_, md):
            return md.get('age', np.inf)

        obs_ids = ['a', 'b', 'c', 'd']
        samp_ids = ['1', '2', '3', '4']
        data = {(0, 0): 1, (0, 1): 2, (0, 2): 3, (0, 3): 4,
                (1, 0): 5, (1, 1): 6, (1, 2): 7, (1, 3): 8,
                (2, 0): 8, (2, 1): 9, (2, 2): 10, (2, 3): 11,
                (3, 0): 12, (3, 1): 13, (3, 2): 14, (3, 3): 15}
        obs_md = [{}, {}, {}, {}]
        samp_md = [{'age': 2, 'foo': 10}, {'age': 4}, {'age': 2, 'bar': 5}, {}]
        t = Table(data, obs_ids, samp_ids, obs_md, samp_md)
        obs_bins, obs_tables = unzip(t.partition(f))

        exp_bins = (2, 4, np.inf)
        exp1_data = {(0, 0): 1, (0, 1): 3, (1, 0): 5, (1, 1): 7, (2, 0): 8,
                     (2, 1): 10, (3, 0): 12, (3, 1): 14}
        exp1_obs_ids = ['a', 'b', 'c', 'd']
        exp1_samp_ids = ['1', '3']
        exp1_obs_md = [{}, {}, {}, {}]
        exp1_samp_md = [{'age': 2, 'foo': 10}, {'age': 2, 'bar': 5}]
        exp1 = Table(exp1_data, exp1_obs_ids, exp1_samp_ids, exp1_obs_md,
                     exp1_samp_md)
        exp2_data = {(0, 0): 2, (1, 0): 6, (2, 0): 9, (3, 0): 13}
        exp2_obs_ids = ['a', 'b', 'c', 'd']
        exp2_samp_ids = ['2']
        exp2_obs_md = [{}, {}, {}, {}]
        exp2_samp_md = [{'age': 4}]
        exp2 = Table(exp2_data, exp2_obs_ids, exp2_samp_ids, exp2_obs_md,
                     exp2_samp_md)
        exp3_data = {(0, 0): 4, (1, 0): 8, (2, 0): 11, (3, 0): 15}
        exp3_obs_ids = ['a', 'b', 'c', 'd']
        exp3_samp_ids = ['4']
        exp3_obs_md = [{}, {}, {}, {}]
        exp3_samp_md = [{}]
        exp3 = Table(exp3_data, exp3_obs_ids, exp3_samp_ids, exp3_obs_md,
                     exp3_samp_md)
        exp_tables = (exp1, exp2, exp3)

        exp1_idx = obs_bins.index(exp_bins[0])
        exp2_idx = obs_bins.index(exp_bins[1])
        exp3_idx = obs_bins.index(exp_bins[2])
        obs_sort = (obs_bins[exp1_idx], obs_bins[exp2_idx], obs_bins[exp3_idx])
        self.assertEqual(obs_sort, exp_bins)
        obs_sort = (obs_tables[exp1_idx], obs_tables[exp2_idx],
                    obs_tables[exp3_idx])

        self.assertEqual(obs_sort, exp_tables)

        # We should get the same table type back.
        exp_types = (Table, Table, Table)
        obs_sort = (type(obs_tables[exp1_idx]), type(obs_tables[exp2_idx]),
                    type(obs_tables[exp3_idx]))
        self.assertEqual(obs_sort, exp_types)

        # Test passing a different constructor. We should get the same data
        # equality, but different table types.
        obs_bins, obs_tables = unzip(t.partition(f))

        obs_sort = (obs_bins[exp1_idx], obs_bins[exp2_idx], obs_bins[exp3_idx])
        self.assertEqual(obs_sort, exp_bins)
        obs_sort = (obs_tables[exp1_idx], obs_tables[exp2_idx],
                    obs_tables[exp3_idx])
        self.assertEqual(obs_sort, exp_tables)
        exp_types = (Table, Table, Table)
        obs_sort = (type(obs_tables[exp1_idx]), type(obs_tables[exp2_idx]),
                    type(obs_tables[exp3_idx]))
        self.assertEqual(obs_sort, exp_types)

    def test_bin_observations_by_metadata(self):
        """Yield tables binned by observation metadata"""
        def make_level_f(level):
            def f(id_, metadata):
                return metadata['taxonomy'][:level]
            return f

        func_king = make_level_f(1)
        func_phy = make_level_f(2)

        obs_ids = ['a', 'b', 'c']
        samp_ids = [1, 2, 3]
        data = {(0, 0): 1, (0, 1): 2, (0, 2): 3,
                (1, 0): 4, (1, 1): 5, (1, 2): 6,
                (2, 0): 7, (2, 1): 8, (2, 2): 9}
        obs_md = [{"taxonomy": ['k__a', 'p__b', 'c__c']},
                  {"taxonomy": ['k__a', 'p__b', 'c__d']},
                  {"taxonomy": ['k__a', 'p__c', 'c__e']}]
        t = Table(data, obs_ids, samp_ids, observation_metadata=obs_md)

        exp_king_obs_ids = ['a', 'b', 'c']
        exp_king_samp_ids = [1, 2, 3]
        exp_king_obs_md = [{"taxonomy": ['k__a', 'p__b', 'c__c']},
                           {"taxonomy": ['k__a', 'p__b', 'c__d']},
                           {"taxonomy": ['k__a', 'p__c', 'c__e']}]
        exp_king = Table(data, exp_king_obs_ids, exp_king_samp_ids,
                         observation_metadata=exp_king_obs_md)
        obs_bins, obs_king = unzip(t.partition(func_king, axis='observation'))

        self.assertEqual(obs_king, [exp_king])
        self.assertEqual(obs_bins, [tuple(['k__a'])])
        self.assertEqual(type(obs_king[0]), type(exp_king))

        obs_bins, obs_king = unzip(t.partition(func_king, axis='observation'))
        self.assertEqual(obs_king, [exp_king])
        self.assertEqual(obs_bins, [tuple(['k__a'])])
        self.assertEqual(type(obs_king[0]), Table)

        exp_phy1_obs_ids = ['a', 'b']
        exp_phy1_samp_ids = [1, 2, 3]
        exp_phy1_data = np.array([[1, 2, 3], [4, 5, 6]])
        exp_phy1_data = {(0, 0): 1, (0, 1): 2, (0, 2): 3,
                         (1, 0): 4, (1, 1): 5, (1, 2): 6}
        exp_phy1_obs_md = [{"taxonomy": ['k__a', 'p__b', 'c__c']},
                           {"taxonomy": ['k__a', 'p__b', 'c__d']}]
        exp_phy1 = Table(exp_phy1_data, exp_phy1_obs_ids, exp_phy1_samp_ids,
                         observation_metadata=exp_phy1_obs_md)
        exp_phy2_obs_ids = ['c']
        exp_phy2_samp_ids = [1, 2, 3]
        exp_phy2_data = {(0, 0): 7, (0, 1): 8, (0, 2): 9}
        exp_phy2_obs_md = [{"taxonomy": ['k__a', 'p__c', 'c__e']}]
        exp_phy2 = Table(exp_phy2_data, exp_phy2_obs_ids, exp_phy2_samp_ids,
                         observation_metadata=exp_phy2_obs_md)
        obs_bins, obs_phy = unzip(t.partition(func_phy, axis='observation'))
        self.assertIn(obs_phy[0], [exp_phy1, exp_phy2])
        self.assertIn(obs_phy[1], [exp_phy1, exp_phy2])
        self.assertIn(obs_bins[0], [('k__a', 'p__b'), ('k__a', 'p__c')])
        self.assertIn(obs_bins[1], [('k__a', 'p__b'), ('k__a', 'p__c')])

    def test_get_table_density(self):
        """Test correctly computes density of table."""
        # Perfectly dense tables.
        npt.assert_almost_equal(self.st1.get_table_density(), 1.0)
        npt.assert_almost_equal(self.st3.get_table_density(), 1.0)
        npt.assert_almost_equal(self.st_rich.get_table_density(), 1.0)

        # Empty table (no dimensions).
        npt.assert_almost_equal(self.empty_st.get_table_density(), 0.0)

        # Tables with some zeros.
        npt.assert_almost_equal(self.st5.get_table_density(), 0.5)

        # Tables with all zeros (with dimensions).
        npt.assert_almost_equal(self.st6.get_table_density(), 0.0)

        # Tables with some zeros explicitly defined.
        npt.assert_almost_equal(self.st7.get_table_density(), 0.75)


class SupportTests2(TestCase):

    def test_coo_arrays_to_sparse(self):
        """convert (values, (row, col)) to scipy"""
        n_rows, n_cols = 3, 4
        exp_d = lil_matrix((n_rows, n_cols))
        exp_d[(0, 0)] = 10
        exp_d[(1, 3)] = 5
        exp_d[(2, 1)] = 2
        exp_d = exp_d.tocoo()
        exp = lil_matrix((n_rows, n_cols))
        exp[(0, 0)] = 10
        exp[(1, 3)] = 5
        exp[(2, 1)] = 2
        data = (np.array([5.0, 2.0, 10.0]),
                (np.array([1, 2, 0]),
                 np.array([3, 1, 0])))
        obs = coo_arrays_to_sparse(data, shape=(n_rows, n_cols))
        self.assertEqual((obs != exp).sum(), 0)

    def test_list_list_to_sparse(self):
        """convert [[row,col,value], ...] to scipy"""
        input = [[0, 0, 1], [1, 1, 5.0], [0, 2, 6]]
        exp = lil_matrix((2, 3))
        exp[(0, 0)] = 1.0
        exp[(1, 1)] = 5.0
        exp[(0, 2)] = 6
        obs = list_list_to_sparse(input)
        self.assertEqual((obs != exp).sum(), 0)

    def test_nparray_to_sparse(self):
        """Convert nparray to sparse"""
        input = np.array([[1, 2, 3, 4], [-1, 6, 7, 8], [9, 10, 11, 12]])
        exp = lil_matrix((3, 4))
        exp[(0, 0)] = 1
        exp[(0, 1)] = 2
        exp[(0, 2)] = 3
        exp[(0, 3)] = 4
        exp[(1, 0)] = -1
        exp[(1, 1)] = 6
        exp[(1, 2)] = 7
        exp[(1, 3)] = 8
        exp[(2, 0)] = 9
        exp[(2, 1)] = 10
        exp[(2, 2)] = 11
        exp[(2, 3)] = 12
        obs = nparray_to_sparse(input)
        self.assertEqual((obs != exp).sum(), 0)

    def test_list_dict_to_sparse(self):
        """Take a list of dicts and condense down to a single dict"""
        input = [{(0, 0): 10, (0, 1): 2}, {(1, 2): 15}, {(0, 3): 7}]
        exp = lil_matrix((3, 4))
        exp[(0, 0)] = 10
        exp[(0, 1)] = 2
        exp[(1, 2)] = 15
        exp[(2, 3)] = 7
        obs = list_dict_to_sparse(input)
        self.assertEqual((obs != exp).sum(), 0)

    def test_dict_to_sparse(self):
        """Take a dict and convert to sparse"""
        input = {(0, 1): 5, (1, 0): 2, (2, 1): 6}
        exp = lil_matrix((3, 2))
        exp[(0, 1)] = 5
        exp[(1, 0)] = 2
        exp[(2, 1)] = 6
        obs = dict_to_sparse(input)
        self.assertEqual((obs != exp).sum(), 0)

    def test_to_sparse(self):
        """Convert to expected sparse types"""
        vals = {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8}
        obs = Table._to_sparse(vals)
        exp = lil_matrix((2, 2))
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(1, 0)] = 7
        exp[(1, 1)] = 8
        self.assertEqual((obs != exp).sum(), 0)

        input = {(0, 1): 5, (10, 8): -1.23}
        input_transpose = {(1, 0): 5, (8, 10): -1.23}

        exp = lil_matrix((11, 9))
        exp[(0, 1)] = 5
        exp[(10, 8)] = -1.23
        obs = Table._to_sparse(input)
        self.assertEqual((obs != exp).sum(), 0)

        # test transpose
        exp = lil_matrix((9, 11))
        exp[(1, 0)] = 5
        exp[(8, 10)] = -1.23
        obs = Table._to_sparse(input_transpose)
        self.assertEqual((obs != exp).sum(), 0)

        # passing a list of dicts, transpose
        exp = lil_matrix((3, 2))
        exp[(0, 0)] = 5.0
        exp[(1, 0)] = 6.0
        exp[(2, 0)] = 7.0
        exp[(0, 1)] = 8.0
        exp[(1, 1)] = 9.0
        exp[(2, 1)] = 10.0
        obs = Table._to_sparse([{(0, 0): 5, (1, 0): 6, (2, 0): 7},
                                {(0, 1): 8, (1, 1): 9, (2, 1): 10}])
        self.assertEqual((obs != exp).sum(), 0)

        # passing a list of lil_matrix
        exp = lil_matrix((2, 3))
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(0, 2)] = 7
        exp[(1, 0)] = 8
        exp[(1, 1)] = 9
        exp[(1, 2)] = 10
        row1 = lil_matrix((1, 3))
        row1[(0, 0)] = 5
        row1[(0, 1)] = 6
        row1[(0, 2)] = 7
        row2 = lil_matrix((1, 3))
        row2[(0, 0)] = 8
        row2[(0, 1)] = 9
        row2[(0, 2)] = 10
        obs = Table._to_sparse([row1, row2])
        self.assertEqual((obs != exp).sum(), 0)

        # test empty set
        exp = lil_matrix((0, 0))
        obs = Table._to_sparse([])
        self.assertEqual((obs != exp).sum(), 0)

    def test_list_nparray_to_sparse(self):
        """lists of nparrays to sparse"""
        ins = [np.array([0, 2, 1, 0]), np.array([1, 0, 0, 1])]
        exp = lil_matrix((2, 4))
        exp[(0, 1)] = 2
        exp[(0, 2)] = 1
        exp[(1, 0)] = 1
        exp[(1, 3)] = 1
        obs = list_nparray_to_sparse(ins)
        self.assertEqual((obs != exp).sum(), 0)

    def test_list_sparse_to_sparse(self):
        """list of lil_matrix to sparse"""
        ins = [lil_matrix((1, 4)), lil_matrix((1, 4))]
        ins[0][0, 0] = 5
        ins[0][0, 1] = 10
        ins[1][0, 2] = 1
        ins[1][0, 3] = 2
        exp = lil_matrix((2, 4))
        exp[0, 0] = 5
        exp[0, 1] = 10
        exp[1, 2] = 1
        exp[1, 3] = 2
        obs = list_sparse_to_sparse(ins)
        self.assertEqual((obs != exp).sum(), 0)


legacy_otu_table1 = """# some comment goes here
#OTU id\tFing\tKey\tNA\tConsensus Lineage
0\t19111\t44536\t42 \tBacteria; Actinobacteria; Actinobacteridae; Propioniba\
cterineae; Propionibacterium

1\t1216\t3500\t6\tBacteria; Firmicutes; Alicyclobacillaceae; Bacilli; La\
ctobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7\t1803\t1184\t2\tBacteria; Actinobacteria; Actinobacteridae; Gordoniace\
ae; Corynebacteriaceae
3\t1722\t4903\t17\tBacteria; Firmicutes; Alicyclobacillaceae; Bacilli; St\
aphylococcaceae
4\t589\t2074\t34\tBacteria; Cyanobacteria; Chloroplasts; vectors
"""
legacy_otu_table_bad_metadata = """# some comment goes here
#OTU id\tFing\tKey\tNA\tConsensus Lineage
0\t19111\t44536\t42 \t
1\t1216\t3500\t6\tBacteria; Firmicutes; Alicyclobacillaceae; Bacilli; La\
ctobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7\t1803\t1184\t2\tBacteria; Actinobacteria; Actinobacteridae; Gordoniace\
ae; Corynebacteriaceae
3\t1722\t4903\t17\tBacteria; Firmicutes; Alicyclobacillaceae; Bacilli; St\
aphylococcaceae
4\t589\t2074\t34\tBacteria; Cyanobacteria; Chloroplasts; vectors
"""
extract_tsv_bug = """#OTU ID	s1	s2	taxonomy
1	123	32\t
2	315	3	k__test;p__test
3	0	22	k__test;p__test"""
otu_table1 = """# Some comment
#OTU ID\tFing\tKey\tNA\tConsensus Lineage
0\t19111\t44536\t42\tBacteria; Actinobacteria; Actinobacteridae; \
Propionibacterineae; Propionibacterium
# some other comment
1\t1216\t3500\t6\tBacteria; Firmicutes; Alicyclobacillaceae; Bacilli; \
Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7\t1803\t1184\t2\tBacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; \
Corynebacteriaceae
# comments
#    everywhere!
3\t1722\t4903\t17\tBacteria; Firmicutes; Alicyclobacillaceae; \
Bacilli; Staphylococcaceae
4\t589\t2074\t34\tBacteria; Cyanobacteria; Chloroplasts; vectors
"""

OBS_META_TYPES = {'sc_separated': lambda x: [e.strip() for e in x.split(';')],
                  'naive': lambda x: x
                  }
OBS_META_TYPES['taxonomy'] = OBS_META_TYPES['sc_separated']

if __name__ == '__main__':
    main()
