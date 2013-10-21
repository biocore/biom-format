#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from numpy import array, float64
from biom.unit_test import TestCase, main
from biom.backends.scipysparse import (ScipySparseMat, to_scipy,
                                       list_nparray_to_scipy,
                                       list_list_to_scipy, list_scipy_to_scipy,
                                       nparray_to_scipy, dict_to_scipy,
                                       list_dict_to_scipy)

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

class ScipySparseMatTests(TestCase):
    def setUp(self):
        # 1 0 2
        # 3 0 4
        self.mat1 = ScipySparseMat(2,3,data=array([[1,0,2],[3,0,4]]))

        # 0x0 case
        self.null = ScipySparseMat(0,0)

        # 0 0
        # 0 0
        self.empty = ScipySparseMat(2,2)

        # 1 0 3
        self.row_vec = ScipySparseMat(1,3)
        self.row_vec[0,0] = 1
        self.row_vec[0,2] = 3

        # 1
        # 0
        # 3
        self.col_vec = ScipySparseMat(3,1)
        self.col_vec[0,0] = 1
        self.col_vec[2,0] = 3

        # 1x1
        self.single_ele = ScipySparseMat(1,1)
        self.single_ele[0,0] = 42

        # Explicit zeros.
        self.explicit_zeros = ScipySparseMat(2,3,
                data=([1,2,3,0,4],([0,0,1,1,1],[0,2,0,1,2])))

    def test_convertVectorToDense(self):
        """Properly converts ScipySparseMat vectors to dense numpy repr."""
        exp = array([1, 0, 3])
        obs = ScipySparseMat.convertVectorToDense(self.row_vec)
        self.assertEqual(obs, exp)

        exp = array([1, 0, 3])
        obs = ScipySparseMat.convertVectorToDense(self.col_vec)
        self.assertEqual(obs, exp)

        exp = array([42])
        obs = ScipySparseMat.convertVectorToDense(self.single_ele)
        self.assertEqual(obs, exp)

    def test_is_empty(self):
        """Differentiate empty (0x0) matrix from non-empty matrix."""
        self.assertTrue(self.null.is_empty)
        self.assertFalse(self.mat1.is_empty)

    def test_shape(self):
        """What kind of shape are you in?"""
        self.assertEqual(self.null.shape, (0,0))
        self.assertEqual(self.mat1.shape, (2,3))
        self.assertEqual(self.empty.shape, (2,2))
        self.assertEqual(self.row_vec.shape, (1,3))
        self.assertEqual(self.col_vec.shape, (3,1))
        self.assertEqual(self.single_ele.shape, (1,1))

    def test_dtype(self):
        """What's your type?"""
        self.assertEqual(self.null.dtype, None)
        self.assertEqual(self.empty.dtype, float)
        self.assertEqual(self.row_vec.dtype, float)

    def test_fmt(self):
        """What format are you in?"""
        self.assertEqual(self.null.fmt, None)
        self.assertEqual(self.empty.fmt, 'coo')
        self.assertEqual(self.mat1.fmt, 'csr')

    def test_size(self):
        """What is your NNZ?"""
        self.assertEqual(self.null.size, 0)
        self.assertEqual(self.empty.size, 0)
        self.assertEqual(self.single_ele.size, 1)
        self.assertEqual(self.mat1.size, 4)
        self.assertEqual(self.explicit_zeros.size, 4)


if __name__ == '__main__':
    main()
