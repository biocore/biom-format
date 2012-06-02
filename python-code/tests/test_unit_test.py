#!/usr/bin/env python

import numpy
from numpy import testing, zeros, array
from biom.unit_test import TestCase, main
from sys import exc_info

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Sandra Smit",
                            "Zongzhi Liu", "Micah Hamady", "Daniel McDonald"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.0.0b"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class TestCaseTests(TestCase):
    """Tests for extension of the built-in unittest framework.

    For each test, includes an example of success and failure.
    """
    unequal_pairs = [
                    (1, 0),
                    ([], ()),
                    (None, 0),
                    ('', ' '),
                    (1, '1'),
                    (0, '0'),
                    ('', None),
                    (array([1,2,3]),array([1,2,4])),
                    (array([[1,2],[3,4]]), array([[1.0,2.0],[3.0,4.1]])),
                    (array([1]), array([1,2])),
                    (zeros(0), array([1])),
                    (array([1,1,1]), array([1])),
                    (array([[1,1],[1,1]]), array([1,1,1,1])),
                    (zeros(0), None),
                    (zeros(3), zeros(5)),
                    (zeros(0), ''),
                ]

    equal_pairs = [
                (1, 1),
                (0, 0),
                (5, 5L),
                (5, 5.0),
                (0, 0.0),
                ('', ''),
                (' ', ' '),
                ('a', 'a'),
                (None, None),
                ([0, 1], [0.0, 1.0]),
                (array([1,2,3]), array([1,2,3])),
                (array([[1,2],[3,4]]), array([[1.0,2.0],[3.0,4.0]])),
                (zeros(0), []),
                (zeros(0), zeros(0)),
                (array([]), zeros(0)),
                (zeros(3), zeros(3)),
                (array([0,0,0]), zeros(3)),
                (array([]), []),
            ]

    small = 1e-7
    big = 1e-5

    within_1e6_abs_pairs = [
                (1, 1 + small),
                (1 + small, 1),
                (1, 1 - small),
                (1 - small, 1),
                (100000, 100000 - small),
                (-100000, -100000 - small),
                (-1, -1 + small),
                (-1, -1 - small),
                (0, small),
                (0, -small),
                (array([1,2]), array([1,2+small])),
                (array([[1,2],[3,4]]), array([[1,2+small],[3,4]]))
                ]

    within_1e6_rel_pairs = [
                (1, 1 + 1 * small),
                (1 + 1 * small, 1),
                (1, 1 - 1 * small),
                (1 - 1 * small, 1),
                (100000, 100000 - 100000 * small),
                (-100000, -100000 - 100000 * small),
                (-1, -1 + -1 * small),
                (-1, -1 - -1 * small),
                (array([1,2]), array([1+small,2])),
                (array([[1000,1000],[1000,1000]]), \
                    array([[1000+1000*small, 1000], [1000,1000]])),
            ]

    outside_1e6_abs_pairs = [
                (1, 1 + big),
                (1 + big, 1),
                (1, 1 - big),
                (1 - big, 1),
                (100000, 100000 - big),
                (-100000, -100000 - big),
                (-1, -1 + big),
                (-1, -1 - big),
                (0, big),
                (0, -big),
                (1e7, 1e7 + 1),
                (array([1,1]), array([1,1+big])),
                (array([[1,1],[1,1]]), array([[1,1+big],[1,1]])),
                ]

    outside_1e6_rel_pairs = [
                (1, 1 + 1 * big),
                (1 + 1 * big, 1),
                (1, 1 - 1 * big),
                (1 - 1 * big, 1),
                (100000, 100000 - 100000 * big),
                (-100000, -100000 - 100000 * big),
                (-1, -1 + -1 * big),
                (-1, -1 - -1 * big),
                (1e-30, 1e-30 + small),
                (0, small),
                (1e5, 1e5 + 1),
                (array([1,1]), array([1,1+1*big])),
            ]

    def test_assertEqual_None(self):
        """assertEqual should not raise exception with two copies of None"""
        try:
            self.assertEqual(None, None)
        except:
            raise AssertionError, \
            "unit_test.assertEqual failed on input %s and %s" \
            % (`first`, `second`)

    def test_assertEqual_numbers(self):
        """assertEqual should not raise exception with integer and float zero"""
        try:
            self.assertEqual(0, 0.0)
        except:
            raise AssertionError, \
            "unit_test.assertEqual failed on input %s and %s" \
            % (`first`, `second`)

    def test_assertEqual_unequal(self):
        """assertEqual should raise exception when values differ"""
        for first, second in self.unequal_pairs:
            try:
                self.assertEqual(first, second)
            except:
                message = str(exc_info()[1])
                self.assertEqual(message,
                'Got %s, but expected %s' \
                % (`first`, `second`))
            else:
                raise AssertionError, \
                "unit_test.assertEqual failed on input %s and %s" \
                % (`first`, `second`)

    def test_assertEqual_equal(self):
        """assertEqual should not raise exception when values test equal"""
        for first, second in self.equal_pairs:
            try:
                self.assertEqual(first, second)
            except:
                raise AssertionError, \
                "unit_test.assertEqual failed on input %s and %s" \
                % (`first`, `second`)

    def test_assertEqual_nested_array(self):
        self.assertEqual([[1,0], [0,1]],
                [array([1,0]), array([0,1])])

    def test_assertEqual_shape_mismatch(self):
        """assertEqual should raise when obs and exp shapes mismatch"""
        obs = [1,2,3]
        exp = [1,2,3,4]
        self.assertRaises(AssertionError, self.assertEqual, obs, exp)

    def test_assertFloatEqualAbs_equal(self):
        """assertFloatEqualAbs should not raise exception when values within eps"""
        for first, second in self.within_1e6_abs_pairs:
            try:
                self.assertFloatEqualAbs(first, second, eps=1e-6)
            except:
                raise AssertionError, \
                "unit_test.assertFloatEqualAbs failed on input %s and %s" \
                % (`first`, `second`)

    def test_assertFloatEqualAbs_threshold(self):
        """assertFloatEqualAbs should raise exception when eps is very small"""
        for first, second in self.within_1e6_abs_pairs:
            try:
                self.assertFloatEqualAbs(first, second, 1e-30)
            except:
                message = str(exc_info()[1])
                diff = first - second
                self.assertEqual(message,
                'Got %s, but expected %s (diff was %s)' \
                % (`first`, `second`, `diff`))
            else:
                raise AssertionError, \
                "unit_test.assertFloatEqualAbs failed on input %s and %s" \
                % (`first`, `second`)


    def test_assertFloatEqualAbs_unequal(self):
        """assertFloatEqualAbs should raise exception when values differ by >eps"""
        for first, second in self.outside_1e6_abs_pairs:
            try:
                self.assertFloatEqualAbs(first, second)
            except:
                message = str(exc_info()[1])
                diff = first - second
                self.assertEqual(message,
                'Got %s, but expected %s (diff was %s)' \
                % (`first`, `second`, `diff`))
            else:
                raise AssertionError, \
                "unit_test.assertFloatEqualAbs failed on input %s and %s" \
                % (`first`, `second`)

    def test_assertFloatEqualAbs_shape_mismatch(self):
        """assertFloatEqualAbs should raise when obs and exp shapes mismatch"""
        obs = [1,2,3]
        exp = [1,2,3,4]
        self.assertRaises(AssertionError, self.assertFloatEqualAbs, obs, exp)

    def test_assertFloatEqualRel_equal(self):
        """assertFloatEqualRel should not raise exception when values within eps"""
        for first, second in self.within_1e6_rel_pairs:
            try:
                self.assertFloatEqualRel(first, second)
            except:
                raise AssertionError, \
                "unit_test.assertFloatEqualRel failed on input %s and %s" \
                % (`first`, `second`)

    def test_assertFloatEqualRel_unequal(self):
        """assertFloatEqualRel should raise exception when eps is very small"""
        for first, second in self.within_1e6_rel_pairs:
            try:
                self.assertFloatEqualRel(first, second, 1e-30)
            except:
                message = str(exc_info()[1])
                diff = first - second
                self.assertEqual(message,
                'Got %s, but expected %s (diff was %s)' \
                % (`first`, `second`, `diff`))
            else:
                raise AssertionError, \
                "unit_test.assertFloatEqualRel failed on input %s and %s" \
                % (`first`, `second`)


    def test_assertFloatEqualRel_unequal(self):
        """assertFloatEqualRel should raise exception when values differ by >eps"""
        for first, second in self.outside_1e6_rel_pairs:
            try:
                self.assertFloatEqualRel(first, second)
            except:
                message = str(exc_info()[1])
                diff = first - second
                self.assertEqual(message,
                'Got %s, but expected %s (diff was %s)' \
                % (`first`, `second`, `diff`))
            else:
                raise AssertionError, \
                "unit_test.assertFloatEqualRel failed on input %s and %s" \
                % (`first`, `second`)

    def test_assertFloatEqualRel_shape_mismatch(self):
        """assertFloatEqualRel should raise when obs and exp shapes mismatch"""
        obs = [1,2,3]
        exp = [1,2,3,4]
        self.assertRaises(AssertionError, self.assertFloatEqualRel, obs, exp)

    def test_assertFloatEqualList_equal(self):
        """assertFloatEqual should work on two lists of similar values"""
        originals = [0, 1, -1, 10, -10, 100, -100]
        modified = [i + 1e-7 for i in originals]
        try:
            self.assertFloatEqual(originals, modified)
            self.assertFloatEqual([], [])   #test empty lists as well
        except:
            raise AssertionError, \
            "unit_test.assertFloatEqual failed on lists of similar values"

    def test_assertFloatEqual_shape_mismatch(self):
        """assertFloatEqual should raise when obs and exp shapes mismatch"""
        obs = [1,2,3]
        exp = [1,2,3,4]
        self.assertRaises(AssertionError, self.assertFloatEqual, obs, exp)

    def test_assertFloatEqualList_unequal(self):
        """assertFloatEqual should fail on two lists of dissimilar values"""
        originals = [0, 1, -1, 10, -10, 100, -100]
        modified = [i + 1e-5 for i in originals]
        try:
            self.assertFloatEqual(originals, modified)
        except:
            pass
        else:
            raise AssertionError, \
            "unit_test.assertFloatEqual failed on lists of dissimilar values"

    def test_assertFloatEqual_mixed(self):
        """assertFloatEqual should work on equal lists of mixed types."""
        first = [i[0] for i in self.equal_pairs]
        second = [i[1] for i in self.equal_pairs]
        self.assertFloatEqual(first, second)

    def test_assertFloatEqualAbs_mixed(self):
        first = [i[0] for i in self.equal_pairs]
        second = [i[1] for i in self.equal_pairs]
        """assertFloatEqualAbs should work on equal lists of mixed types."""
        self.assertFloatEqualAbs(first, second)

    def test_assertFloatEqualRel_mixed(self):
        first = [i[0] for i in self.equal_pairs]
        second = [i[1] for i in self.equal_pairs]
        """assertFloatEqualRel should work on equal lists of mixed types."""
        self.assertFloatEqualRel(first, second)

    def test_assertFloatEqual_mixed_unequal(self):
        """assertFloatEqual should work on unequal lists of mixed types."""
        first = [i[0] for i in self.unequal_pairs]
        second = [i[1] for i in self.unequal_pairs]
        self.assertRaises(AssertionError, \
            self.assertFloatEqual, first, second)

    def test_assertFloatEqualAbs_mixed(self):
        """assertFloatEqualAbs should work on lists of mixed types."""
        first = [i[0] for i in self.unequal_pairs]
        second = [i[1] for i in self.unequal_pairs]
        self.assertRaises(AssertionError, \
            self.assertFloatEqualAbs, first, second)

    def test_assertFloatEqualRel_mixed(self):
        """assertFloatEqualRel should work on lists of mixed types."""
        first = [i[0] for i in self.unequal_pairs]
        second = [i[1] for i in self.unequal_pairs]
        self.assertRaises(AssertionError, \
            self.assertFloatEqualRel, first, second)

if __name__ == '__main__':
    main()
