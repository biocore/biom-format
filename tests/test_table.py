#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from tempfile import mktemp
from unittest import TestCase, main

import h5py
import numpy.testing as npt
from numpy import where, zeros, array

from biom.util import unzip
from biom.table import (TableException, Table, UnknownID,
                        prefer_self, index_list, dict_to_nparray,
                        list_dict_to_nparray, table_factory,
                        list_list_to_nparray, to_sparse,
                        nparray_to_sparseobj, list_nparray_to_sparseobj,
                        SparseObj, get_zerod_matrix)
from biom.parse import parse_biom_table_hdf5

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Justin Kuczynski",
               "Greg Caporaso", "Jose Clemente", "Adam Robbins-Pianka"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class SupportTests(TestCase):
    def test_get_zerod_matrix(self):
        """returns a zerod matrix"""
        foo = array([[1, 2, 3], [4, 5, 6]])
        exp = zeros((2, 3))
        obs = get_zerod_matrix(foo)
        npt.assert_equal(obs, exp)

        foo = SparseObj(2, 3)
        foo[1, 2] = 3
        exp = SparseObj(2, 3)
        obs = get_zerod_matrix(foo)
        self.assertEqual(obs, exp)

    def test_table_factory_sparseobj_nparray(self):
        """beat the table_factory sparsely to death"""
        # nparray test
        samp_ids = ['1', '2', '3', '4']
        obs_ids = ['a', 'b', 'c']
        nparray = array([[1, 2, 3, 4], [-1, 6, 7, 8], [9, 10, 11, 12]])
        data = nparray_to_sparseobj(
            array([[1, 2, 3, 4], [-1, 6, 7, 8], [9, 10, 11, 12]]))
        exp = Table(data, samp_ids, obs_ids)
        obs = table_factory(nparray, samp_ids, obs_ids)
        self.assertEqual(obs, exp)

    def test_table_factory_sparse_list_nparray(self):
        """beat the table_factory sparsely to death"""
        # list of nparray test
        samp_ids = ['1', '2', '3', '4']
        obs_ids = ['a', 'b', 'c']
        list_np = [array([1, 2, 3, 4]), array([5, 6, 7, 8]),
                   array([9, 10, 11, 12])]
        data = list_nparray_to_sparseobj(list_np)
        exp = Table(data, samp_ids, obs_ids)
        obs = table_factory(list_np, samp_ids, obs_ids)
        self.assertEqual(obs, exp)

    def test_table_factory_sparse_dict(self):
        """beat the table_factory sparsely to death"""
        # dict test
        samp_ids = range(24)
        obs_ids = range(101)
        dict_input = {(0, 0): 1, (0, 10): 5, (100, 23): -3}
        d_input = zeros((101, 24), dtype=float)
        d_input[0, 0] = 1
        d_input[0, 10] = 5
        d_input[100, 23] = -3
        data = nparray_to_sparseobj(d_input)
        exp = Table(data, samp_ids, obs_ids)
        obs = table_factory(dict_input, samp_ids, obs_ids)
        self.assertEqual(obs, exp)

    def test_table_factory_sparse_list_dict(self):
        """beat the table_factory sparsely to death"""
        # list of dict test
        samp_ids = range(11)
        obs_ids = range(3)
        ld_input = zeros((3, 11), dtype=float)
        ld_input[0, 5] = 10
        ld_input[0, 10] = 2
        ld_input[1, 1] = 15
        ld_input[2, 3] = 7
        data = nparray_to_sparseobj(ld_input)
        exp = Table(data, samp_ids, obs_ids)
        list_dict = [{(0, 5): 10, (10, 10): 2}, {(0, 1): 15}, {(0, 3): 7}]
        obs = table_factory(list_dict, samp_ids, obs_ids)
        self.assertEqual(obs, exp)

    def test_table_factory_sparse_list_list(self):
        """beat the table_factory sparsely to death"""
        # list list test
        samp_ids = range(3)
        obs_ids = range(2)
        exp_data = SparseObj(2, 3)
        exp_data[0, 1] = 5
        exp_data[1, 2] = 10
        exp = Table(exp_data, samp_ids, obs_ids)
        input_ = [[0, 1, 5], [1, 2, 10]]
        obs = table_factory(input_, samp_ids, obs_ids)
        self.assertEqual(obs, exp)

    def test_table_exception(self):
        """Make sure a TableException can be raised"""
        def f():
            raise TableException
        self.assertRaises(TableException, f)

    def test_list_list_to_nparray(self):
        """Convert [[value, value, ... value], ...] to nparray"""
        input_ = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], [7, 6, 5, 4, 3]]
        exp = array(
            [[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], [7, 6, 5, 4, 3]], dtype=float)
        obs = list_list_to_nparray(input_)
        npt.assert_equal(obs, exp)

    def test_dict_to_nparray(self):
        """Take a dict -> array"""
        input_ = {(0, 0): 1, (0, 10): 5, (100, 23): -3}
        exp = zeros((101, 24), dtype=float)
        exp[0, 0] = 1
        exp[0, 10] = 5
        exp[100, 23] = -3
        obs = dict_to_nparray(input_)
        npt.assert_equal(obs, exp)

    def test_list_dict_to_nparray(self):
        """List of dict -> nparray"""
        input_ = [{(0, 5): 10, (10, 10): 2}, {(0, 1): 15}, {(0, 3): 7}]
        exp = zeros((3, 11), dtype=float)
        exp[0, 5] = 10
        exp[0, 10] = 2
        exp[1, 1] = 15
        exp[2, 3] = 7
        obs = list_dict_to_nparray(input_)
        npt.assert_equal(obs, exp)

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
            to_sparse(array([[5, 6], [7, 8]])), [1, 2], [3, 4])
        self.vals = {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8}
        self.st1 = Table(to_sparse(self.vals),
                         ['a', 'b'], ['1', '2'])
        self.st2 = Table(to_sparse(self.vals),
                         ['a', 'b'], ['1', '2'])
        self.vals3 = to_sparse({(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4})
        self.vals4 = to_sparse({(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4})
        self.st3 = Table(self.vals3, ['b', 'c'], ['2', '3'])
        self.st4 = Table(self.vals4, ['c', 'd'], ['3', '4'])
        self._to_dict_f = lambda x: sorted(x.items())
        self.st_rich = Table(to_sparse(self.vals),
                             ['a', 'b'], ['1', '2'],
                             [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                             [{'taxonomy': ['k__a', 'p__b']},
                              {'taxonomy': ['k__a', 'p__c']}])

        self.empty_st = Table(to_sparse([]), [], [])

        self.vals5 = to_sparse({(0, 1): 2, (1, 1): 4})
        self.st5 = Table(self.vals5, ['a', 'b'], ['5', '6'])

        self.vals6 = to_sparse({(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0})
        self.st6 = Table(self.vals6, ['a', 'b'], ['5', '6'])

        self.vals7 = to_sparse({(0, 0): 5, (0, 1): 7, (1, 0): 8, (1, 1): 0})
        self.st7 = Table(self.vals7, ['a', 'b'], ['5', '6'])

        self.single_sample_st = Table(
            to_sparse(array([[2.0], [0.0], [1.0]])), ['S1'],
            ['O1', 'O2', 'O3'])
        self.single_obs_st = Table(to_sparse(array([[2.0, 0.0, 1.0]])),
                                   ['S1', 'S2', 'S3'], ['O1'])

        self.to_remove = []

    def tearDown(self):
        if self.to_remove:
            for f in self.to_remove:
                os.remove(f)

    def test_format_hdf5(self):
        """Write a file"""
        fname = mktemp()
        self.to_remove.append(fname)
        h5 = h5py.File(fname, 'w')
        self.st_rich.format_hdf5(h5, 'tests')
        h5.close()

        h5 = h5py.File(fname, 'r')
        self.assertIn('observation', h5)
        self.assertIn('sample', h5)
        self.assertEqual(sorted(h5.attrs.keys()), sorted(['id', 'type',
                                                          'format-url',
                                                          'format-version',
                                                          'generated-by',
                                                          'creation-date',
                                                          'shape', 'nnz']))
        obs = parse_biom_table_hdf5(h5)
        self.assertEqual(obs, self.st_rich)

    def test_get_sample_index(self):
        """returns the sample index"""
        self.assertEqual(0, self.simple_derived.get_sample_index(1))
        self.assertEqual(1, self.simple_derived.get_sample_index(2))
        self.assertRaises(UnknownID, self.simple_derived.get_sample_index, 3)

    def test_get_observation_index(self):
        """returns the observation index"""
        self.assertEqual(0, self.simple_derived.get_observation_index(3))
        self.assertEqual(1, self.simple_derived.get_observation_index(4))
        self.assertRaises(
            UnknownID,
            self.simple_derived.get_observation_index,
            5)

    def test_index_ids(self):
        """Index the all the ids!!!"""
        exp_samp = {1: 0, 2: 1}
        exp_obs = {3: 0, 4: 1}
        self.assertEqual(self.simple_derived._sample_index, exp_samp)
        self.assertEqual(self.simple_derived._obs_index, exp_obs)

    def test_sample_exists(self):
        """Verify samples exist!"""
        self.assertTrue(self.simple_derived.sample_exists(1))
        self.assertTrue(self.simple_derived.sample_exists(2))
        self.assertFalse(self.simple_derived.sample_exists(3))

    def test_observation_exists(self):
        """Verify observation exist!"""
        self.assertTrue(self.simple_derived.observation_exists(3))
        self.assertTrue(self.simple_derived.observation_exists(4))
        self.assertFalse(self.simple_derived.observation_exists(2))

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
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        Table(d, samp_ids, obs_ids, samp_md, obs_md)
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
        Table(d, samp_ids, obs_ids, samp_md, obs_md)

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
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md, obs_md)

        self.assertEqual(t.sample_metadata[0]['non existent key'], None)
        self.assertEqual(t.sample_metadata[1]['non existent key'], None)
        self.assertEqual(t.sample_metadata[2]['non existent key'], None)
        self.assertEqual(t.sample_metadata[3]['non existent key'], None)
        self.assertEqual(t.observation_metadata[0]['non existent key'], None)
        self.assertEqual(t.observation_metadata[1]['non existent key'], None)
        self.assertEqual(t.observation_metadata[2]['non existent key'], None)

    def test_add_observation_metadata_w_existing_metadata(self):
        """ add_observationMetadata functions with existing metadata """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 9}, {'a': 8}, {'a': 7}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 0}, {'e': 0}, {'f': 0}, {'g': 0}]
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md, obs_md)
        self.assertEqual(t.observation_metadata[0]['a'], 9)
        self.assertEqual(t.observation_metadata[1]['a'], 8)
        self.assertEqual(t.observation_metadata[2]['a'], 7)
        obs_md = {1: {'taxonomy': ['A', 'B']},
                  2: {'taxonomy': ['B', 'C']},
                  3: {'taxonomy': ['E', 'D', 'F']},
                  4: {'taxonomy': ['this', 'is', 'ignored']}}
        t.add_observation_metadata(obs_md)
        self.assertEqual(t.observation_metadata[0]['a'], 9)
        self.assertEqual(t.observation_metadata[1]['a'], 8)
        self.assertEqual(t.observation_metadata[2]['a'], 7)
        self.assertEqual(t.observation_metadata[0]['taxonomy'], ['A', 'B'])
        self.assertEqual(t.observation_metadata[1]['taxonomy'], ['B', 'C'])
        self.assertEqual(t.observation_metadata[2]['taxonomy'],
                                                  ['E', 'D', 'F'])

    def test_add_observation_metadata_one_entry(self):
        """ add_observationMetadata functions with single md entry """
        obs_ids = [1, 2, 3]
        obs_md = {1: {'taxonomy': ['A', 'B']},
                  2: {'taxonomy': ['B', 'C']},
                  3: {'taxonomy': ['E', 'D', 'F']}}
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 0}, {'e': 0}, {'f': 0}, {'g': 0}]
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md, obs_md=None)
        t.add_observation_metadata(obs_md)
        self.assertEqual(t.observation_metadata[0]['taxonomy'], ['A', 'B'])
        self.assertEqual(t.observation_metadata[1]['taxonomy'], ['B', 'C'])
        self.assertEqual(t.observation_metadata[2]['taxonomy'],
                                                  ['E', 'D', 'F'])

    def test_add_observation_metadata_two_entries(self):
        """ add_observationMetadata functions with more than one md entry """
        obs_ids = [1, 2, 3]
        obs_md = {1: {'taxonomy': ['A', 'B'], 'other': 'h1'},
                  2: {'taxonomy': ['B', 'C'], 'other': 'h2'},
                  3: {'taxonomy': ['E', 'D', 'F'], 'other': 'h3'}}
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'d': 0}, {'e': 0}, {'f': 0}, {'g': 0}]
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md, obs_md=None)
        t.add_observation_metadata(obs_md)
        self.assertEqual(t.observation_metadata[0]['taxonomy'], ['A', 'B'])
        self.assertEqual(t.observation_metadata[1]['taxonomy'], ['B', 'C'])
        self.assertEqual(t.observation_metadata[2]['taxonomy'],
                         ['E', 'D', 'F'])
        self.assertEqual(t.observation_metadata[0]['other'], 'h1')
        self.assertEqual(t.observation_metadata[1]['other'], 'h2')
        self.assertEqual(t.observation_metadata[2]['other'], 'h3')

    def test_add_sample_metadata_one_w_existing_metadata(self):
        """ add_sample_metadata functions with existing metadata """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = [{'Treatment': 'Control'},
                   {'Treatment': 'Fasting'},
                   {'Treatment': 'Fasting'},
                   {'Treatment': 'Control'}]
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md, obs_md=obs_md)
        self.assertEqual(t.sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t.sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[3]['Treatment'], 'Control')

        samp_md = {4: {'barcode': 'TTTT'},
                   6: {'barcode': 'AAAA'},
                   5: {'barcode': 'GGGG'},
                   7: {'barcode': 'CCCC'},
                   10: {'ignore': 'me'}}
        t.add_sample_metadata(samp_md)
        self.assertEqual(t.sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t.sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[3]['Treatment'], 'Control')
        self.assertEqual(t.sample_metadata[0]['barcode'], 'TTTT')
        self.assertEqual(t.sample_metadata[1]['barcode'], 'GGGG')
        self.assertEqual(t.sample_metadata[2]['barcode'], 'AAAA')
        self.assertEqual(t.sample_metadata[3]['barcode'], 'CCCC')

    def test_add_sample_metadata_one_entry(self):
        """ add_sample_metadata functions with single md entry """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = {4: {'Treatment': 'Control'},
                   5: {'Treatment': 'Fasting'},
                   6: {'Treatment': 'Fasting'},
                   7: {'Treatment': 'Control'}}
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md=None, obs_md=obs_md)
        t.add_sample_metadata(samp_md)
        self.assertEqual(t.sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t.sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[3]['Treatment'], 'Control')

    def test_add_sample_metadata_two_entries(self):
        """ add_sample_metadata functions with more than one md entry """
        obs_ids = [1, 2, 3]
        obs_md = [{'a': 0}, {'b': 0}, {'c': 0}]
        samp_ids = [4, 5, 6, 7]
        samp_md = {4: {'Treatment': 'Control', 'D': ['A', 'A']},
                   5: {'Treatment': 'Fasting', 'D': ['A', 'B']},
                   6: {'Treatment': 'Fasting', 'D': ['A', 'C']},
                   7: {'Treatment': 'Control', 'D': ['A', 'D']}}
        d = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        t = Table(d, samp_ids, obs_ids, samp_md=None, obs_md=obs_md)
        t.add_sample_metadata(samp_md)
        self.assertEqual(t.sample_metadata[0]['Treatment'], 'Control')
        self.assertEqual(t.sample_metadata[1]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[2]['Treatment'], 'Fasting')
        self.assertEqual(t.sample_metadata[3]['Treatment'], 'Control')
        self.assertEqual(t.sample_metadata[0]['D'], ['A', 'A'])
        self.assertEqual(t.sample_metadata[1]['D'], ['A', 'B'])
        self.assertEqual(t.sample_metadata[2]['D'], ['A', 'C'])
        self.assertEqual(t.sample_metadata[3]['D'], ['A', 'D'])

    def test_get_value_by_ids(self):
        """Return the value located in the matrix by the ids"""
        t1 = Table(array([[5, 6], [7, 8]]), [1, 2], [3, 4])
        t2 = Table(array([[5, 6], [7, 8]]), ['a', 'b'], ['c', 'd'])

        self.assertEqual(5, t1.get_value_by_ids(3, 1))
        self.assertEqual(6, t1.get_value_by_ids(3, 2))
        self.assertEqual(7, t1.get_value_by_ids(4, 1))
        self.assertEqual(8, t1.get_value_by_ids(4, 2))
        self.assertEqual(5, t2.get_value_by_ids('c', 'a'))
        self.assertEqual(6, t2.get_value_by_ids('c', 'b'))
        self.assertEqual(7, t2.get_value_by_ids('d', 'a'))
        self.assertEqual(8, t2.get_value_by_ids('d', 'b'))

        self.assertRaises(UnknownID, t1.get_value_by_ids, 'a', 1)
        self.assertRaises(UnknownID, t2.get_value_by_ids, 0, 0)

    def test_getitem(self):
        """getitem should work as expeceted"""
        self.assertEqual(self.simple_derived[0, 0], 5)
        self.assertEqual(self.simple_derived[1, 0], 7)
        self.assertEqual(self.simple_derived[0, 1], 6)
        self.assertEqual(self.simple_derived[1, 1], 8)
        self.assertRaises(IndexError, self.simple_derived.__getitem__, [1, 2])

    def test_is_empty(self):
        """returns true if empty"""
        self.assertTrue(Table(array([]), [], []).is_empty())
        self.assertFalse(self.simple_derived.is_empty())

    def test_immutability(self):
        """Test Table object immutability."""
        # Try to set members to something completely different.
        self.assertRaises(TypeError, self.st1.__setattr__, 'sample_ids',
                          ['foo', 'bar'])
        self.assertRaises(TypeError, self.simple_derived.__setattr__,
                          'observation_ids', ['foo', 'bar'])
        self.assertRaises(TypeError, self.st2.__setattr__,
                          'sample_metadata', [{'foo': 'bar'}, {'bar': 'baz'}])
        self.assertRaises(TypeError, self.st1.__setattr__,
                          'observation_metadata', [{'foo': 'bar'},
                                                   {'bar': 'baz'}])


class SparseTableTests(TestCase):

    def setUp(self):
        self.vals = {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8}
        self.st1 = Table(to_sparse(self.vals),
                         ['a', 'b'], ['1', '2'])
        self.st2 = Table(to_sparse(self.vals),
                         ['a', 'b'], ['1', '2'])
        self.vals3 = to_sparse({(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4})
        self.vals4 = to_sparse({(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4})
        self.st3 = Table(self.vals3, ['b', 'c'], ['2', '3'])
        self.st4 = Table(self.vals4, ['c', 'd'], ['3', '4'])
        self._to_dict_f = lambda x: sorted(x.items())
        self.st_rich = Table(to_sparse(self.vals),
                             ['a', 'b'], ['1', '2'],
                             [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                             [{'taxonomy': ['k__a', 'p__b']},
                              {'taxonomy': ['k__a', 'p__c']}])

        self.empty_st = Table(to_sparse([]), [], [])

        self.vals5 = to_sparse({(0, 1): 2, (1, 1): 4})
        self.st5 = Table(self.vals5, ['a', 'b'], ['5', '6'])

        self.vals6 = to_sparse({(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0})
        self.st6 = Table(self.vals6, ['a', 'b'], ['5', '6'])

        self.vals7 = to_sparse({(0, 0): 5, (0, 1): 7, (1, 0): 8, (1, 1): 0})
        self.st7 = Table(self.vals7, ['a', 'b'], ['5', '6'])

        self.single_sample_st = Table(
            to_sparse(array([[2.0], [0.0], [1.0]])), ['S1'],
            ['O1', 'O2', 'O3'])
        self.single_obs_st = Table(to_sparse(array([[2.0, 0.0, 1.0]])),
                                   ['S1', 'S2', 'S3'], ['O1'])

    def test_sum(self):
        """Test of sum!"""
        self.assertEqual(self.st1.sum('whole'), 26)
        npt.assert_equal(self.st1.sum('sample'), array([12, 14]))
        npt.assert_equal(self.st1.sum('observation'), array([11, 15]))

        exp = array([3.0])
        obs = self.single_sample_st.sum('sample')
        self.assertEqual(obs, exp)

        exp = array([3.0])
        obs = self.single_obs_st.sum('observation')
        self.assertEqual(obs, exp)

    def test_reduce(self):
        """Reduce method"""
        f = lambda x, y: x * 2 + y
        npt.assert_equal(self.st1.reduce(f, 'sample'), array([17, 20]))
        npt.assert_equal(self.st1.reduce(f, 'observation'), array([16, 22]))

    def test_transpose(self):
        """Should transpose a sparse table"""
        obs = self.st1.transpose()

        self.assertEqual(obs.sample_ids, self.st1.observation_ids)
        self.assertEqual(obs.observation_ids, self.st1.sample_ids)
        npt.assert_equal(obs.sample_data('1'), self.st1.observation_data('1'))
        npt.assert_equal(obs.sample_data('2'), self.st1.observation_data('2'))
        self.assertEqual(obs.transpose(), self.st1)

        obs = self.st_rich.transpose()

        self.assertEqual(obs.sample_ids, self.st_rich.observation_ids)
        self.assertEqual(obs.observation_ids, self.st_rich.sample_ids)
        self.assertEqual(obs.sample_metadata,
                         self.st_rich.observation_metadata)
        self.assertEqual(obs.observation_metadata,
                         self.st_rich.sample_metadata)
        npt.assert_equal(obs.sample_data('1'),
                         self.st_rich.observation_data('1'))
        npt.assert_equal(obs.sample_data('2'),
                         self.st_rich.observation_data('2'))
        self.assertEqual(obs.transpose(), self.st_rich)

    def test_sort_observation_order(self):
        """sort by observations arbitrary order"""
        vals = {(0, 0): 7, (0, 1): 8, (1, 0): 5, (1, 1): 6}
        exp = Table(to_sparse(vals),
                    ['a', 'b'], ['2', '1'])
        obs = self.st1.sort_observation_order(['2', '1'])
        self.assertEqual(obs, exp)

    def test_sort_sample_order(self):
        """sort by observations arbitrary order"""
        vals = {(0, 0): 6, (0, 1): 5,
                (1, 0): 8, (1, 1): 7}
        exp = Table(to_sparse(vals),
                    ['b', 'a'], ['1', '2'])

        obs = self.st1.sort_sample_order(['b', 'a'])
        self.assertEqual(obs, exp)

    def test_sort_by_sample_id(self):
        """sort by samples by a function"""
        sort_f = sorted
        data_in = nparray_to_sparseobj(
            array([[1, 2, 3, 8], [4, 5, 6, 9], [7, 8, 9, 11]]))
        t = Table(data_in, ['c', 'a', 'b', 'd'], [2, 1, 3])
        exp_data = nparray_to_sparseobj(
            array([[2, 3, 1, 8], [5, 6, 4, 9], [8, 9, 7, 11]]))
        exp = Table(exp_data, ['a', 'b', 'c', 'd'], [2, 1, 3])
        obs = t.sort_by_sample_id(sort_f=sort_f)
        self.assertEqual(obs, exp)

    def test_sort_by_observation_id(self):
        """sort by observation ids by a function"""
        sort_f = sorted
        data_in = nparray_to_sparseobj(
            array([[1, 2, 3, 8], [4, 5, 6, 9], [7, 8, 9, 11]]), float)
        t = Table(data_in, ['c', 'a', 'b', 'd'], [2, 1, 3])
        exp_data = nparray_to_sparseobj(
            array([[4, 5, 6, 9], [1, 2, 3, 8], [7, 8, 9, 11]]), float)
        exp = Table(exp_data, ['c', 'a', 'b', 'd'], [1, 2, 3])
        obs = t.sort_by_observation_id(sort_f=sort_f)
        self.assertEqual(obs, exp)

    def test_eq(self):
        """sparse equality"""
        self.assertTrue(self.st1 == self.st2)
        object.__setattr__(self.st1, 'observation_ids', [1, 2, 3])
        self.assertFalse(self.st1 == self.st2)

        object.__setattr__(self.st1, 'observation_ids',
                           self.st2.observation_ids)
        object.__setattr__(self.st1, '_data',
                           nparray_to_sparseobj(array([[1, 2], [10, 20]])))
        self.assertFalse(self.st1 == self.st2)

    def test_data_equality(self):
        """check equality between tables"""
        self.assertTrue(self.st1._data_equality(self.st2))
        self.assertTrue(self.st1._data_equality(self.st1))
        self.assertFalse(self.st1._data_equality(self.st3))

    def test_nonzero(self):
        """Return a list of nonzero positions"""
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(to_sparse(data), ['a', 'b', 'c', 'd'], ['1', '2', '3'])
        exp = [('1', 'a'), ('1', 'b'), ('1', 'd'), ('2', 'b'), ('2', 'd'),
               ('3', 'a'), ('3', 'b')]
        obs = list(st.nonzero())
        self.assertEqual(obs, exp)

    def test_nonzero_counts(self):
        """Returns nonzero counts over an axis"""
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(to_sparse(data), ['a', 'b', 'c', 'd'], ['1', '2', '3'])

        exp_samp = array([6, 12, 0, 11])
        exp_obs = array([14, 15, 0])
        exp_whole = array([29])

        obs_samp = st.nonzero_counts('sample')
        obs_obs = st.nonzero_counts('observation')
        obs_whole = st.nonzero_counts('whole')

        npt.assert_equal(obs_samp, exp_samp)
        npt.assert_equal(obs_obs, exp_obs)
        npt.assert_equal(obs_whole, exp_whole)

    def test_nonzero_counts_binary(self):
        """Returns nonzero counts over an axis"""
        data = {(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 3,
                (1, 0): 0, (1, 1): 7, (1, 2): 0, (1, 3): 8,
                (2, 0): 1, (2, 1): -1, (2, 2): 0, (2, 3): 0}
        st = Table(to_sparse(data), ['a', 'b', 'c', 'd'], ['1', '2', '3'])

        exp_samp = array([2, 3, 0, 2])
        exp_obs = array([3, 2, 2])
        exp_whole = array([7])

        obs_samp = st.nonzero_counts('sample', binary=True)
        obs_obs = st.nonzero_counts('observation', binary=True)
        obs_whole = st.nonzero_counts('whole', binary=True)

        npt.assert_equal(obs_samp, exp_samp)
        npt.assert_equal(obs_obs, exp_obs)
        npt.assert_equal(obs_whole, exp_whole)

    def test_merge(self):
        """Merge two tables"""
        u = 'union'
        i = 'intersection'

        # test 1
        data = to_sparse({(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16})
        exp = Table(data, ['a', 'b'], ['1', '2'])
        obs = self.st1.merge(self.st1, sample=u, observation=u)
        self.assertEqual(obs, exp)

        # test 2
        data = to_sparse(
            {(0, 0): 5, (0, 1): 6, (0, 2): 0, (1, 0): 7, (1, 1): 9, (1, 2): 2,
             (2, 0): 0, (2, 1): 3, (2, 2): 4})
        exp = Table(data, ['a', 'b', 'c'], ['1', '2', '3'])
        obs = self.st1.merge(self.st3, sample=u, observation=u)
        self.assertEqual(obs, exp)

        # test 3
        data = to_sparse({(0, 0): 5, (0, 1): 6, (0, 2): 0, (0, 3): 0,
                          (1, 0): 7, (1, 1): 8, (1, 2): 0, (1, 3): 0,
                          (2, 0): 0, (2, 1): 0, (2, 2): 1, (2, 3): 2,
                          (3, 0): 0, (3, 1): 0, (3, 2): 3, (3, 3): 4})
        exp = Table(data, ['a', 'b', 'c', 'd'], ['1', '2', '3', '4'])
        obs = self.st1.merge(self.st4, sample=u, observation=u)
        self.assertEqual(obs, exp)

        # test 4
        data = to_sparse({(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16})
        exp = Table(data, ['a', 'b'], ['1', '2'])
        obs = self.st1.merge(self.st1, sample=i, observation=i)
        self.assertEqual(obs, exp)

        # test 5
        exp = Table(to_sparse({(0, 0): 9}), ['b'], ['2'])
        obs = self.st1.merge(self.st3, sample=i, observation=i)
        self.assertEqual(obs, exp)

        # test 6
        self.assertRaises(TableException, self.st1.merge, self.st4, i, i)

        # test 7
        data = to_sparse({(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16})
        exp = Table(data, ['a', 'b'], ['1', '2'])
        obs = self.st1.merge(self.st1, sample=i, observation=u)
        self.assertEqual(obs, exp)

        # test 8
        data = to_sparse({(0, 0): 6, (1, 0): 9, (2, 0): 3})
        exp = Table(data, ['b'], ['1', '2', '3'])
        obs = self.st1.merge(self.st3, sample=i, observation=u)
        self.assertEqual(obs, exp)

        # test 9
        self.assertRaises(TableException, self.st1.merge, self.st4, i, u)

        # test 10
        data = to_sparse({(0, 0): 10, (0, 1): 12, (1, 0): 14, (1, 1): 16})
        exp = Table(data, ['a', 'b'], ['1', '2'])
        obs = self.st1.merge(self.st1, sample=u, observation=i)
        self.assertEqual(obs, exp)

        # test 11
        data = to_sparse({(0, 0): 7, (0, 1): 9, (0, 2): 2})
        exp = Table(data, ['a', 'b', 'c'], ['2'])
        obs = self.st1.merge(self.st3, sample=u, observation=i)
        self.assertEqual(obs, exp)

        # test 12
        self.assertRaises(TableException, self.st1.merge, self.st4, u, i)

    def test_sample_data(self):
        """tested in derived class"""
        exp = array([5, 7])
        obs = self.st1.sample_data('a')
        npt.assert_equal(obs, exp)
        self.assertRaises(UnknownID, self.st1.sample_data, 'asdasd')

    def test_observation_data(self):
        """tested in derived class"""
        exp = array([5, 6])
        obs = self.st1.observation_data('1')
        npt.assert_equal(obs, exp)
        self.assertRaises(UnknownID, self.st1.observation_data, 'asdsad')

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

    def test_conv_to_np(self):
        """Should convert a self styled vector to numpy type"""
        input_row = SparseObj(1, 3)
        input_row[(0, 0)] = 10
        exp = array([10.0, 0, 0])
        obs = self.st1._conv_to_np(input_row)
        npt.assert_equal(obs, exp)

        input_col = SparseObj(3, 1)
        input_col[(0, 0)] = 12
        exp = array([12.0, 0, 0])
        obs = self.st1._conv_to_np(input_col)
        npt.assert_equal(obs, exp)

        # 1x1
        input_vec = SparseObj(1, 1)
        input_vec[(0, 0)] = 42
        exp = array([42.0])
        obs = self.st1._conv_to_np(input_vec)
        npt.assert_equal(obs, exp)

    def test_conv_to_self_type(self):
        """Should convert other to SparseObj type"""
        exp = SparseObj(2, 2)
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(1, 0)] = 7
        exp[(1, 1)] = 8
        obs = self.st1._conv_to_self_type(self.vals)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        exp = SparseObj(2, 2)
        exp[(0, 0)] = 5
        exp[(0, 1)] = 7
        exp[(1, 0)] = 6
        exp[(1, 1)] = 8
        obs = self.st1._conv_to_self_type(self.vals, transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a single vector
        exp = SparseObj(1, 3)
        exp[(0, 0)] = 2
        exp[(0, 1)] = 0
        exp[(0, 2)] = 3
        obs = self.st1._conv_to_self_type(array([2, 0, 3]))
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of dicts
        exp = SparseObj(2, 3)
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(0, 2)] = 7
        exp[(1, 0)] = 8
        exp[(1, 1)] = 9
        exp[(1, 2)] = 10
        obs = self.st1._conv_to_self_type([{(0, 0): 5, (0, 1): 6, (0, 2): 7},
                                           {(1, 0): 8, (1, 1): 9, (1, 2): 10}])
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

    def test_iter(self):
        """Should iterate over samples"""
        exp = [(array([5, 7]), 'a', None), (array([6, 8]), 'b', None)]
        obs = list(self.st1)
        npt.assert_equal(obs, exp)

    def test_iter_obs(self):
        """Iterate over observations of sparse matrix"""
        r1 = SparseObj(1, 2)
        r2 = SparseObj(1, 2)
        r1[(0, 0)] = 5
        r1[(0, 1)] = 6
        r2[(0, 0)] = 7
        r2[(0, 1)] = 8

        exp = map(self._to_dict_f, [r1, r2])
        obs = map(self._to_dict_f, self.st1._iter_obs())

        self.assertEqual(obs, exp)

    def test_iter_samp(self):
        """Iterate over samples of sparse matrix"""
        c1 = SparseObj(1, 2)
        c2 = SparseObj(1, 2)
        c1[(0, 0)] = 5
        c1[(0, 1)] = 7
        c2[(0, 0)] = 6
        c2[(0, 1)] = 8
        exp = map(self._to_dict_f, [c1, c2])
        obs = map(self._to_dict_f, self.st1._iter_samp())

        self.assertEqual(obs, exp)

    def test_iter_samples(self):
        """Iterates samples"""
        gen = self.st1.iter_samples()
        exp = [(array([5, 7]), 'a', None), (array([6, 8]), 'b', None)]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter_samples()
        exp = [(array([5, 7]), 'a', {'barcode': 'aatt'}),
               (array([6, 8]), 'b', {'barcode': 'ttgg'})]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        # [[1,2,3],[1,0,2]] isn't yielding column 2 correctly
        vals = {(0, 0): 5, (0, 1): 6, (1, 1): 8}
        st = Table(to_sparse(vals), ['a', 'b'], ['1', '2'])
        gen = st.iter_samples()
        exp = [(array([5, 0]), 'a', None), (array([6, 8]), 'b', None)]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_observations(self):
        """Iterates observations"""
        gen = self.st1.iter_observations()
        exp = [(array([5, 6]), '1', None), (array([7, 8]), '2', None)]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter_observations()
        exp = [(array([5, 6]), '1', {'taxonomy': ['k__a', 'p__b']}),
               (array([7, 8]), '2', {'taxonomy': ['k__a', 'p__c']})]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_sample_data(self):
        """Iterates data by samples"""
        gen = self.st1.iter_sample_data()
        exp = [array([5, 7]), array([6, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter_sample_data()
        exp = [array([5, 7]), array([6, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        # [[1,2,3],[1,0,2]] isn't yielding column 2 correctly
        vals = {(0, 0): 5, (0, 1): 6, (1, 1): 8}
        st = Table(to_sparse(vals), ['a', 'b'], ['1', '2'])
        gen = st.iter_sample_data()
        exp = [array([5, 0]), array([6, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_sample_data_single_obs(self):
        """Iterates data by samples with a single observation."""
        exp = [array([2.0]), array([0.0]), array([1.0])]
        obs = list(self.single_obs_st.iter_sample_data())
        # We test this way to make sure the observed value is a single element
        # array instead of a numpy scalar.
        for o, e in zip(obs, exp):
            self.assertEqual(o, e)

    def test_iter_observation_data(self):
        """Iterates data by observations"""
        gen = self.st1.iter_observation_data()
        exp = [array([5, 6]), array([7, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

        gen = self.st_rich.iter_observation_data()
        exp = [array([5, 6]), array([7, 8])]
        obs = list(gen)
        npt.assert_equal(obs, exp)

    def test_iter_observation_data_single_sample(self):
        """Iterates data by observations from a single sample."""
        exp = [array([2.0]), array([0.0]), array([1.0])]
        obs = list(self.single_sample_st.iter_observation_data())
        for o, e in zip(obs, exp):
            self.assertEqual(o, e)

    def test_filter_samples(self):
        """Filters samples by arbitrary function"""
        f_value = lambda v, id_, md: (v <= 5).any()
        f_id = lambda v, id_, md: id_ == 'a'
        f_md = lambda v, id_, md: md['barcode'] == 'ttgg'

        val_sd = to_sparse({(0, 0): 5, (1, 0): 7})
        exp_value = Table(val_sd, ['a'], ['1', '2'],
                          [{'barcode': 'aatt'}],
                          [{'taxonomy': ['k__a', 'p__b']},
                           {'taxonomy': ['k__a', 'p__c']}])
        id_sd = to_sparse({(0, 0): 5, (1, 0): 7})
        exp_id = Table(id_sd, ['a'], ['1', '2'],
                       [{'barcode': 'aatt'}],
                       [{'taxonomy': ['k__a', 'p__b']},
                        {'taxonomy': ['k__a', 'p__c']}])
        md_sd = to_sparse({(0, 0): 6, (1, 0): 8})
        exp_md = Table(md_sd, ['b'], ['1', '2'],
                       [{'barcode': 'ttgg'}], [{'taxonomy': ['k__a', 'p__b']},
                                               {'taxonomy': ['k__a', 'p__c']}])

        obs_value = self.st_rich.filter_samples(f_value)
        obs_id = self.st_rich.filter_samples(f_id)
        obs_md = self.st_rich.filter_samples(f_md)

        self.assertEqual(obs_value, exp_value)
        self.assertEqual(obs_id, exp_id)
        self.assertEqual(obs_md, exp_md)

        inv_sd = to_sparse({(0, 0): 6, (1, 0): 8})
        exp_inv = Table(inv_sd, ['b'], ['1', '2'],
                        [{'barcode': 'ttgg'}],
                        [{'taxonomy': ['k__a', 'p__b']},
                         {'taxonomy': ['k__a', 'p__c']}])
        obs_inv = self.st_rich.filter_samples(f_value, invert=True)
        self.assertEqual(obs_inv, exp_inv)
        self.assertRaises(TableException, self.st_rich.filter_samples,
                          lambda x, y, z: False)

    def test_filter_observations(self):
        """Filters observations by arbitrary function"""
        f_value = lambda v, id_, md: (v <= 5).any()
        f_id = lambda v, id_, md: id_ == '1'
        f_md = lambda v, id_, md: md['taxonomy'][1] == 'p__c'

        val_sd = to_sparse({(0, 0): 5, (0, 1): 6})
        exp_value = Table(val_sd, ['a', 'b'], ['1'],
                          [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                          [{'taxonomy': ['k__a', 'p__b']}])
        id_sd = to_sparse({(0, 0): 5, (0, 1): 6})
        exp_id = Table(id_sd, ['a', 'b'], ['1'],
                       [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                       [{'taxonomy': ['k__a', 'p__b']}])
        md_sd = to_sparse({(0, 0): 7, (0, 1): 8})
        exp_md = Table(md_sd, ['a', 'b'], ['2'],
                       [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                       [{'taxonomy': ['k__a', 'p__c']}])

        obs_value = self.st_rich.filter_observations(f_value)
        obs_id = self.st_rich.filter_observations(f_id)
        obs_md = self.st_rich.filter_observations(f_md)

        self.assertEqual(obs_value, exp_value)
        self.assertEqual(obs_id, exp_id)
        self.assertEqual(obs_md, exp_md)

        inv_sd = to_sparse({(0, 0): 7, (0, 1): 8})
        exp_inv = Table(inv_sd, ['a', 'b'], ['2'],
                        [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                        [{'taxonomy': ['k__a', 'p__c']}])
        obs_inv = self.st_rich.filter_observations(f_value, invert=True)
        self.assertEqual(obs_inv, exp_inv)
        self.assertRaises(TableException, self.st_rich.filter_observations,
                          lambda x, y, z: False)

    def test_transform_observations(self):
        """Transform observations by arbitrary function"""
        def transform_f(v, id, md):
            return where(v >= 7, 1, 0)
        sp_sd = to_sparse({(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 1})
        exp = Table(sp_sd, ['a', 'b'], ['1', '2'])
        obs = self.st1.transform_observations(transform_f)
        self.assertEqual(obs, exp)

    def test_transform_samples(self):
        """Transform samples by arbitrary function"""
        def transform_f(v, id, md):
            return where(v >= 6, 1, 0)
        sp_sd = to_sparse({(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 1})
        exp = Table(sp_sd, ['a', 'b'], ['1', '2'])
        obs = self.st1.transform_samples(transform_f)
        self.assertEqual(obs, exp)

    def test_norm_observation_by_sample(self):
        """normalize observations by sample"""
        data = to_sparse({(0, 0): 2, (0, 1): 0, (1, 0): 6, (1, 1): 1})
        data_exp = to_sparse(
            {(0, 0): 0.25, (0, 1): 0.0, (1, 0): 0.75, (1, 1): 1.0})
        st = Table(data, ['a', 'b'], ['1', '2'])
        exp = Table(data_exp, ['a', 'b'], ['1', '2'])
        obs = st.norm_observation_by_sample()
        self.assertEqual(obs, exp)

    def test_norm_observation_by_metadata(self):
        """normalize observations by sample"""
        data = to_sparse({(0, 0): 6, (0, 1): 0, (1, 0): 6, (1, 1): 1})
        data_exp = to_sparse(
            {(0, 0): 2., (0, 1): 0.0, (1, 0): 3.0, (1, 1): 0.5})
        st = Table(data, ['a', 'b'], ['1', '2'],
                   [{}, {}], [{'CopyNumber': 3}, {'CopyNumber': 2}])
        exp = Table(data_exp, ['a', 'b'], ['1', '2'],
                    [{}, {}], [{'CopyNumber': 3}, {'CopyNumber': 2}])
        obs = st.norm_observation_by_metadata('CopyNumber')
        self.assertEqual(obs, exp)

    def test_norm_sample_by_observation(self):
        """normalize sample by observation"""
        data = to_sparse({(0, 0): 0, (0, 1): 2, (1, 0): 2, (1, 1): 6})
        data_exp = to_sparse(
            {(0, 0): 0.0, (0, 1): 1.0, (1, 0): 0.25, (1, 1): 0.75})
        st = Table(data, ['a', 'b'], ['1', '2'])
        exp = Table(data_exp, ['a', 'b'], ['1', '2'])
        obs = st.norm_sample_by_observation()
        self.assertEqual(obs, exp)

    def test_bin_samples_by_metadata(self):
        """Yield tables binned by sample metadata"""
        f = lambda x: x['age']
        obs_ids = ['a', 'b', 'c', 'd']
        samp_ids = ['1', '2', '3', '4']
        data = {(0, 0): 1, (0, 1): 2, (0, 2): 3, (0, 3): 4,
                (1, 0): 5, (1, 1): 6, (1, 2): 7, (1, 3): 8,
                (2, 0): 8, (2, 1): 9, (2, 2): 10, (2, 3): 11,
                (3, 0): 12, (3, 1): 13, (3, 2): 14, (3, 3): 15}
        obs_md = [{}, {}, {}, {}]
        samp_md = [{'age': 2, 'foo': 10}, {'age': 4}, {'age': 2, 'bar': 5}, {}]
        t = Table(to_sparse(data), samp_ids, obs_ids, samp_md, obs_md)
        obs_bins, obs_tables = unzip(t.bin_samples_by_metadata(f))

        exp_bins = (2, 4, None)
        exp1_data = to_sparse(
            {(0, 0): 1, (0, 1): 3, (1, 0): 5, (1, 1): 7, (2, 0): 8,
             (2, 1): 10, (3, 0): 12, (3, 1): 14})
        exp1_obs_ids = ['a', 'b', 'c', 'd']
        exp1_samp_ids = ['1', '3']
        exp1_obs_md = [{}, {}, {}, {}]
        exp1_samp_md = [{'age': 2, 'foo': 10}, {'age': 2, 'bar': 5}]
        exp1 = Table(exp1_data, exp1_samp_ids, exp1_obs_ids, exp1_samp_md,
                     exp1_obs_md)
        exp2_data = to_sparse({(0, 0): 2, (1, 0): 6, (2, 0): 9, (3, 0): 13})
        exp2_obs_ids = ['a', 'b', 'c', 'd']
        exp2_samp_ids = ['2']
        exp2_obs_md = [{}, {}, {}, {}]
        exp2_samp_md = [{'age': 4}]
        exp2 = Table(exp2_data, exp2_samp_ids, exp2_obs_ids, exp2_samp_md,
                     exp2_obs_md)
        exp3_data = to_sparse({(0, 0): 4, (1, 0): 8, (2, 0): 11, (3, 0): 15})
        exp3_obs_ids = ['a', 'b', 'c', 'd']
        exp3_samp_ids = ['4']
        exp3_obs_md = [{}, {}, {}, {}]
        exp3_samp_md = [{'age': None}]
        exp3 = Table(exp3_data, exp3_samp_ids, exp3_obs_ids, exp3_samp_md,
                     exp3_obs_md)
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
        obs_bins, obs_tables = unzip(t.bin_samples_by_metadata(f))

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
            def f(metadata):
                return metadata['taxonomy'][:level]
            return f

        func_king = make_level_f(1)
        func_phy = make_level_f(2)

        obs_ids = ['a', 'b', 'c']
        samp_ids = [1, 2, 3]
        data = to_sparse({(0, 0): 1, (0, 1): 2, (0, 2): 3,
                          (1, 0): 4, (1, 1): 5, (1, 2): 6,
                          (2, 0): 7, (2, 1): 8, (2, 2): 9})
        obs_md = [{"taxonomy": ['k__a', 'p__b', 'c__c']},
                  {"taxonomy": ['k__a', 'p__b', 'c__d']},
                  {"taxonomy": ['k__a', 'p__c', 'c__e']}]
        t = Table(data, samp_ids, obs_ids, observation_metadata=obs_md)

        exp_king_obs_ids = ['a', 'b', 'c']
        exp_king_samp_ids = [1, 2, 3]
        exp_king_obs_md = [{"taxonomy": ['k__a', 'p__b', 'c__c']},
                           {"taxonomy": ['k__a', 'p__b', 'c__d']},
                           {"taxonomy": ['k__a', 'p__c', 'c__e']}]
        exp_king = Table(data, exp_king_samp_ids, exp_king_obs_ids,
                         observation_metadata=exp_king_obs_md)
        obs_bins, obs_king = unzip(t.bin_observations_by_metadata(func_king))

        self.assertEqual(obs_king, [exp_king])
        self.assertEqual(obs_bins, [tuple(['k__a'])])
        self.assertEqual(type(obs_king[0]), type(exp_king))

        obs_bins, obs_king = unzip(t.bin_observations_by_metadata(func_king))
        self.assertEqual(obs_king, [exp_king])
        self.assertEqual(obs_bins, [tuple(['k__a'])])
        self.assertEqual(type(obs_king[0]), Table)

        exp_phy1_obs_ids = ['a', 'b']
        exp_phy1_samp_ids = [1, 2, 3]
        exp_phy1_data = array([[1, 2, 3], [4, 5, 6]])
        exp_phy1_data = to_sparse({(0, 0): 1, (0, 1): 2, (0, 2): 3,
                                   (1, 0): 4, (1, 1): 5, (1, 2): 6})
        exp_phy1_obs_md = [{"taxonomy": ['k__a', 'p__b', 'c__c']},
                           {"taxonomy": ['k__a', 'p__b', 'c__d']}]
        exp_phy1 = Table(exp_phy1_data, exp_phy1_samp_ids,
                         exp_phy1_obs_ids,
                         observation_metadata=exp_phy1_obs_md)
        exp_phy2_obs_ids = ['c']
        exp_phy2_samp_ids = [1, 2, 3]
        exp_phy2_data = to_sparse({(0, 0): 7, (0, 1): 8, (0, 2): 9})
        exp_phy2_obs_md = [{"taxonomy": ['k__a', 'p__c', 'c__e']}]
        exp_phy2 = Table(exp_phy2_data, exp_phy2_samp_ids, exp_phy2_obs_ids,
                         observation_metadata=exp_phy2_obs_md)
        obs_bins, obs_phy = unzip(t.bin_observations_by_metadata(func_phy))
        self.assertEqual(obs_phy, [exp_phy1, exp_phy2])
        self.assertEqual(obs_bins, [('k__a', 'p__b'), ('k__a', 'p__c')])

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


class SparseOTUTableTests(TestCase):

    def setUp(self):
        self.vals = {(0, 0): 5, (1, 0): 7, (1, 1): 8}
        self.sot_min = Table(
            to_sparse(self.vals, dtype=int), ['a', 'b'], ['1', '2'])
        self.sot_rich = Table(to_sparse(self.vals, dtype=int),
                              ['a', 'b'], ['1', '2'],
                              [{'barcode': 'aatt'}, {'barcode': 'ttgg'}],
                              [{'taxonomy': ['k__a', 'p__b']},
                               {'taxonomy': ['k__a', 'p__c']}])
        self.float_table = Table(to_sparse({(0, 1): 2.5, (0, 2): 3.4,
                                            (1, 0): 9.3, (1, 1): 10.23,
                                            (1, 2): 2.2}),
                                 ['a', 'b', 'c'], ['1', '2'])

    def test_get_biom_format_object_no_generated_by(self):
        """Should raise without a generated_by string"""
        self.assertRaises(
            TableException,
            self.sot_min.get_biom_format_object,
            None)
        self.assertRaises(TableException,
                          self.sot_min.get_biom_format_object, 10)

    def test_get_biom_format_object_minimal(self):
        """Should return a dictionary of the minimal table in Biom format."""
        exp = {'rows': [{'id': '1', 'metadata': None},
                        {'id': '2', 'metadata': None}],
               'format': 'Biological Observation Matrix 1.0.0',
               'data': [[0, 0, 5.0], [1, 0, 7.0], [1, 1, 8.0]],
               'columns': [{'id': 'a', 'metadata': None},
                           {'id': 'b', 'metadata': None}],
               'shape': [2, 2],
               'format_url': __url__,
               'id': None,
               'generated_by': 'foo',
               'matrix_element_type': 'int'}
        obs = self.sot_min.get_biom_format_object('foo')
        del obs['date']
        self.assertEqual(obs, exp)

    def test_get_biom_format_object_rich(self):
        """Should return a dictionary of the rich table in Biom format."""
        exp = {
            'rows': [{'id': '1', 'metadata': {'taxonomy': ['k__a', 'p__b']}},
                     {'id': '2', 'metadata': {'taxonomy': ['k__a', 'p__c']}}],
            'format': 'Biological Observation Matrix 1.0.0',
            'data': [[0, 0, 5.0], [1, 0, 7.0], [1, 1, 8.0]],
            'columns': [{'id': 'a', 'metadata': {'barcode': 'aatt'}},
                        {'id': 'b', 'metadata': {'barcode': 'ttgg'}}],
            'shape': [2, 2],
            'format_url': __url__,
            'id': None,
            'generated_by': 'foo',
            'matrix_element_type': 'int'}
        obs = self.sot_rich.get_biom_format_object('foo')
        del obs['date']
        self.assertEqual(obs, exp)

    def test_get_biom_format_object_float(self):
        """Should return a dictionary of the table with float values."""
        exp = {'rows': [{'id': '1', 'metadata': None},
                        {'id': '2', 'metadata': None}],
               'format': 'Biological Observation Matrix 1.0.0',
               'data': [[0, 1, 2.5],
                        [0, 2, 3.3999999999999999],
                        [1, 0, 9.3000000000000007],
                        [1, 1, 10.23],
                        [1, 2, 2.2000000000000002]],
               'columns': [{'id': 'a', 'metadata': None},
                           {'id': 'b', 'metadata': None},
                           {'id': 'c', 'metadata': None}],
               'shape': [2, 3],
               'format_url': __url__,
               'generated_by': 'foo',
               'id': None,
               'matrix_element_type': 'float'}
        obs = self.float_table.get_biom_format_object('foo')
        del obs['date']
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
