#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from copy import deepcopy

import numpy as np

from biom import example_table, Table
from biom.exception import TableException
from biom.err import (_zz_test_empty, _test_obssize, _test_sampsize,
                      _test_obsdup, _test_sampdup, _test_obsmdsize,
                      _test_sampmdsize, errstate, geterr, seterr, geterrcall,
                      _test_hasnan, _test_hasinf, seterrcall, errcheck,
                      __errprof, IGNORE, RAISE, EMPTY, OBSSIZE, SAMPSIZE, CALL,
                      WARN, OBSDUP, SAMPDUP, OBSMDSIZE, SAMPMDSIZE, HASNAN,
                      HASINF)


runtime_ep = __errprof
runtime_ep_profile = deepcopy(runtime_ep._profile)
runtime_ep_state = runtime_ep._state.copy()
runtime_ep_test = runtime_ep._test.copy()


class ErrModeTests(TestCase):
    def setUp(self):
        self.ex_table = example_table.copy()

    def test_test_empty(self):
        self.assertTrue(_zz_test_empty(Table([], [], [])))
        self.assertFalse(_zz_test_empty(self.ex_table))

    def test_test_obssize(self):
        self.assertFalse(_test_obssize(self.ex_table))
        self.ex_table._observation_ids = self.ex_table._observation_ids[:-1]
        self.assertTrue(_test_obssize(self.ex_table))

    def test_test_sampsize(self):
        self.assertFalse(_test_sampsize(self.ex_table))
        self.ex_table._sample_ids = self.ex_table._sample_ids[:-1]
        self.assertTrue(_test_sampsize(self.ex_table))

    def test_test_obsdup(self):
        self.assertFalse(_test_obsdup(self.ex_table))
        self.ex_table._observation_ids[0] = self.ex_table._observation_ids[1]
        self.assertTrue(_test_obsdup(self.ex_table))

    def test_test_sampdup(self):
        self.assertFalse(_test_sampdup(self.ex_table))
        self.ex_table._sample_ids[0] = self.ex_table._sample_ids[1]
        self.assertTrue(_test_sampdup(self.ex_table))

    def test_test_obsmdsize(self):
        self.assertFalse(_test_obsdup(self.ex_table))
        self.ex_table._observation_metadata = \
            self.ex_table._observation_metadata[:-1]
        self.assertTrue(_test_obsmdsize(self.ex_table))

    def test_test_sampmdsize(self):
        self.assertFalse(_test_sampdup(self.ex_table))
        self.ex_table._sample_metadata = \
            self.ex_table._sample_metadata[:-1]
        self.assertTrue(_test_sampmdsize(self.ex_table))

    def test_test_hasnan(self):
        self.assertFalse(_test_hasnan(self.ex_table))
        self.ex_table._data.data[0] = np.nan
        self.assertTrue(_test_hasnan(self.ex_table))

    def test_test_hasinf(self):
        self.assertFalse(_test_hasinf(self.ex_table))
        self.ex_table._data.data[0] = np.inf
        self.assertTrue(_test_hasinf(self.ex_table))


class ErrorProfileTests(TestCase):
    def setUp(self):
        self.ex_table = example_table.copy()
        self.ep = runtime_ep
        self.ep.state = {'all': 'raise'}

    def tearDown(self):
        self.ep._profile = deepcopy(runtime_ep_profile.copy())
        self.ep._state = runtime_ep_state.copy()
        self.ep._test = runtime_ep_test.copy()

    def test_test(self):
        self.ep.test(self.ex_table)
        self.ep.test(self.ex_table, 'empty')
        self.ep.test(self.ex_table, 'empty', 'obssize')

        self.ex_table._observation_ids = self.ex_table._observation_ids[:-1]
        self.ep.test(self.ex_table, 'empty')
        self.assertTrue(isinstance(self.ep.test(self.ex_table, 'obssize'),
                                   TableException))

    def test_test_evaluation_order(self):
        # issue 813
        tab = Table(np.array([[1, 2], [3, 4]]), ['A', 'B'], ['C', 'D'])
        tab._observation_ids = np.array(['A', 'A'], dtype='object')
        tab._sample_ids = np.array(['B', 'B'], dtype='object')

        self.assertEqual(self.ep.test(tab, 'obsdup', 'sampdup').args[0],
                         'Duplicate observation IDs')
        self.assertEqual(self.ep.test(tab, 'sampdup', 'obsdup').args[0],
                         'Duplicate observation IDs')

    def test_state(self):
        self.ep.state = {'all': IGNORE}
        self.assertEqual(set(self.ep._state.values()), {'ignore'})
        self.ep.state = {'empty': CALL}
        self.assertEqual(set(self.ep._state.values()), {'ignore', CALL})
        self.assertEqual(self.ep.state['empty'], CALL)

        with self.assertRaises(KeyError):
            self.ep.state = {'empty': 'missing'}

        with self.assertRaises(KeyError):
            self.ep.state = {'emptyasdasd': 'ignore'}

    def test_contains(self):
        self.assertTrue('empty' in self.ep)
        self.assertFalse('emptyfoo' in self.ep)

    def test_handle_error(self):
        def callback(foo):
            return 10

        self.ep.setcall('empty', callback)

        self.assertTrue(isinstance(self.ep._handle_error('empty', None),
                        TableException))

        self.ep.state = {'empty': CALL}
        self.assertEqual(self.ep._handle_error('empty', None), 10)

    def test_setcall(self):
        def callback(foo):
            return 10

        self.assertEqual(self.ep._profile['empty'][CALL](None), None)
        self.ep.setcall('empty', callback)
        self.assertEqual(self.ep._profile['empty'][CALL](None), 10)

        with self.assertRaises(KeyError):
            self.ep.setcall('emptyfoo', callback)

    def test_getcall(self):
        def callback(foo):
            return 10
        self.ep.setcall('empty', callback)
        self.assertEqual(self.ep.getcall('empty'), callback)

        with self.assertRaises(KeyError):
            self.ep.getcall('emptyfoo')

    def test_register_unregister(self):
        def cb(x):
            return 123

        def test(x):
            return x == 5

        self.ep.register('foo', 'bar', IGNORE, test, callback=cb)
        self.assertTrue('foo' in self.ep)
        self.ep.state = {'foo': CALL}
        self.assertEqual(self.ep._handle_error('foo', None), 123)

        foo_prof = self.ep._profile['foo'].copy()
        prof, func, state = self.ep.unregister('foo')

        self.assertEqual(func, test)
        self.assertEqual(state, CALL)
        self.assertEqual(prof, foo_prof)

        with self.assertRaises(KeyError):
            self.ep.register('empty', 1, 2, lambda: None)

        with self.assertRaises(KeyError):
            self.ep.register('foo', 'missing', 2, lambda: None)

        with self.assertRaises(KeyError):
            self.ep.unregister('non_existant')


class SupportTests(TestCase):
    def setUp(self):
        self.ex_table = example_table.copy()

    def test_geterr(self):
        state = geterr()
        self.assertEqual(state, runtime_ep._state)
        old = seterr(all=CALL)
        self.assertNotEqual(geterr(), state)
        seterr(**old)

    def test_seterr(self):
        existing = seterr(empty=WARN)
        self.assertEqual(runtime_ep._state['empty'], WARN)
        self.assertNotEqual(runtime_ep._state['empty'], existing)
        seterr(empty=existing['empty'])
        self.assertNotEqual(runtime_ep._state['empty'], WARN)
        self.assertEqual(runtime_ep._state, existing)

    def test_geterrcall(self):
        exp = runtime_ep._profile['sampsize'][CALL]
        obs = geterrcall('sampsize')
        self.assertEqual(obs, exp)

        with self.assertRaises(KeyError):
            geterrcall('asdasd')

    def test_seterrcall(self):
        def foo(x):
            return 10

        seterrcall('sampmdsize', foo)
        obs = geterrcall('sampmdsize')
        self.assertEqual(obs, foo)

    def test_errcheck(self):
        self.assertEqual(errcheck(self.ex_table), None)
        self.ex_table._sample_ids = self.ex_table._sample_ids[:-1]
        with self.assertRaises(TableException):
            errcheck(self.ex_table)

    def test_errstate(self):
        def foo(item):
            return "the callback called"

        table = Table([], [], [])
        seterrcall('empty', foo)
        self.assertNotEqual(geterr()['empty'], CALL)
        with errstate(empty=CALL):
            result = errcheck(table)
        self.assertEqual(result, "the callback called")
        self.assertNotEqual(geterr()['empty'], CALL)

def _what_to_raise(errtype):
    d = {k: IGNORE for k in __errprof._state}
    d[errtype] = RAISE
    return d


class IntegrationTests(TestCase):
    def _check(self, errcond, msg, table_data):
        with self.assertRaisesRegex(TableException, msg):
            with errstate(**_what_to_raise(errcond)):
                Table(*table_data)

    def test_has_duplicate_samples(self):
        data = (np.array([[1, 2, 3], [4, 5, 6]]),
                list('ab'),
                ['S1', 'S1', 'S2'])
        self._check('sampdup', SAMPDUP, data)

    def test_has_duplicate_observations(self):
        data = (np.array([[1, 2, 3], [4, 5, 6]]),
                ['x', 'x'],
                list('abc'))
        self._check('obsdup', OBSDUP, data)

    def test_is_empty(self):
        data = ([], [], [])
        self._check('empty', EMPTY, data)

    def test_observation_size(self):
        data = (np.array([[1, 2, 3], [4, 5, 6]]),
                ['w', 'x', 'y'],
                list('abc'))
        self._check('obssize', OBSSIZE, data)

    def test_sample_size(self):
        data = (np.array([[1, 2, 3], [4, 5, 6]]),
                ['w', 'x'],
                list('ab'))
        self._check('sampsize', SAMPSIZE, data)

    def test_observation_metadata_size(self):
        data = (np.array([[1, 2, 3], [4, 5, 6]]),
                ['x', 'y'],
                list('abc'),
                [{1: 2}, {1: 3}, {1: 4}])
        self._check('obsmdsize', OBSMDSIZE, data)

    def test_sample_metadata_size(self):
        data = (np.array([[1, 2, 3], [4, 5, 6]]),
                ['x', 'y'],
                list('abc'),
                None,
                [{1: 2}, ])
        self._check('sampmdsize', SAMPMDSIZE, data)

    def test_has_nan(self):
        data = (np.array([[1, 2, np.nan], [4, 5, 6]]),
                ['x', 'y'],
                list('abc'))
        self._check('hasnan', HASNAN, data)

    def test_has_inf(self):
        data = (np.array([[1, 2, np.inf], [4, 5, 6]]),
                ['x', 'y'],
                list('abc'))
        self._check('hasinf', HASINF, data)


if __name__ == '__main__':
    main()
