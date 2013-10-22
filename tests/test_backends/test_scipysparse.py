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
        self.assertEqual(self.single_ele.fmt, 'lil')

    def test_size(self):
        """What is your NNZ?"""
        self.assertEqual(self.null.size, 0)
        self.assertEqual(self.empty.size, 0)
        self.assertEqual(self.single_ele.size, 1)
        self.assertEqual(self.mat1.size, 4)
        self.assertEqual(self.explicit_zeros.size, 4)

    def test_convert(self):
        """Test sparse format conversion."""
        self.assertEqual(self.mat1.fmt, 'csr')
        self.mat1.convert('coo')
        self.assertEqual(self.mat1.fmt, 'coo')

        self.assertEqual(self.null.fmt, None)
        self.null.convert('coo')
        self.assertEqual(self.null.fmt, None)

    def test_transpose(self):
        """Test transposition."""
        obs = self.null.T
        self.assertEqual(obs, self.null)
        self.assertFalse(obs is self.null)

        obs = self.single_ele.T
        self.assertEqual(obs, self.single_ele)
        self.assertFalse(obs is self.single_ele)

        exp = ScipySparseMat(3,2,data=array([[1,3],[0,0],[2,4]]))
        obs = self.mat1.T
        self.assertEqual(obs, exp)

    def test_sum(self):
        """Test summing a matrix."""
        self.assertEqual(self.null.sum(), 0)
        self.assertEqual(self.mat1.sum(), 10)
        self.assertEqual(self.mat1.sum(0), array([4,0,6]))
        self.assertEqual(self.mat1.sum(1), array([3,7]))
        self.assertEqual(self.row_vec.sum(1), array([4]))
        self.assertEqual(self.col_vec.sum(0), array([4]))
        with self.assertRaises(ValueError):
            _ = self.mat1.sum(3)

    def test_getRow(self):
        """Test grabbing a row from the matrix."""
        with self.assertRaises(IndexError):
            _ = self.null.getRow(0)

        exp = ScipySparseMat(1,3,data=array([[1,0,2]]))
        obs = self.mat1.getRow(0)
        self.assertEqual(obs, exp)

        self.assertEqual(self.row_vec.getRow(0), self.row_vec)

    def test_getCol(self):
        """Test grabbing a column from the matrix."""
        with self.assertRaises(IndexError):
            _ = self.null.getCol(0)

        exp = ScipySparseMat(2,1,data=array([[1],[3]]))
        obs = self.mat1.getCol(0)
        self.assertEqual(obs, exp)

        self.assertEqual(self.col_vec.getCol(0), self.col_vec)

    def test_items_iteritems(self):
        """Test getting a list of non-zero elements."""
        exp = []
        self.assertEqual(self.null.items(), exp)
        self.assertEqual(list(self.null.iteritems()), exp)

        exp = []
        self.assertEqual(self.empty.items(), exp)
        self.assertEqual(list(self.empty.iteritems()), exp)

        exp = [((0,0),1),((0,2),2),((1,0),3),((1,2),4)]
        self.assertEqual(sorted(self.mat1.items()), exp)
        self.assertEqual(sorted(self.mat1.iteritems()), exp)

    def test_copy(self):
        """Test copying the matrix."""
        copy = self.mat1.copy()
        self.assertEqual(copy, self.mat1)
        self.assertFalse(copy is self.mat1)

        copy[1,1] = 42
        self.assertNotEqual(copy, self.mat1)

    def test_eq(self):
        """Test whether two matrices are equal."""
        self.assertTrue(self.null == ScipySparseMat(0,0))
        self.assertTrue(self.empty == ScipySparseMat(2,2))

        mat2 = ScipySparseMat(2,3,data=array([[1,0,2],[3,0,4]]))
        self.assertTrue(self.mat1 == mat2)

        # Sparse format shouldn't matter.
        mat2.convert('lil')
        self.assertNotEqual(self.mat1.fmt, mat2.fmt)
        self.assertTrue(self.mat1 == mat2)

        # Equality works in both directions.
        self.assertTrue(mat2 == self.mat1)

    def test_ne(self):
        """Test whether two matrices are not equal."""
        # Wrong type.
        self.assertTrue(self.null != array([]))

        # Wrong shape.
        self.assertTrue(self.empty != ScipySparseMat(2,1))

        # Wrong dtype.
        self.assertTrue(self.empty != ScipySparseMat(2,2,dtype=int))

        # Wrong size.
        wrong_size = ScipySparseMat(2,2)
        self.assertTrue(self.empty == wrong_size)
        wrong_size[1,0] = 42
        self.assertTrue(self.empty != wrong_size)

        # Wrong size.
        wrong_data = self.mat1.copy()
        self.assertTrue(self.mat1 == wrong_data)
        wrong_data[0,2] = 42
        self.assertTrue(self.mat1 != wrong_data)
        self.assertTrue(wrong_data != self.mat1)

    def test_str(self):
        """Test getting string representation of the matrix."""
        self.assertEqual(str(self.null), '<0x0 sparse matrix>')

    def test_setitem(self):
        """Test setting an element in the matrix."""
        with self.assertRaises(IndexError):
            self.null[0,0] = 42

        with self.assertRaises(IndexError):
            self.empty[0] = [42,42]

        with self.assertRaises(ValueError):
            self.mat1[0,0] = 0

        with self.assertRaises(ValueError):
            self.mat1[0,0] = 0.0

        # Setting existing zero element doesn't change matrix.
        copy = self.mat1.copy()
        copy[0,1] = 0.0
        self.assertEqual(copy, self.mat1)

        # Setting existing nonzero element to the same thing doesn't change
        # matrix.
        copy = self.mat1.copy()
        copy[0,0] = 1.0
        self.assertEqual(copy, self.mat1)

        # nonzero element -> nonzero element
        copy = self.mat1.copy()
        copy[0,0] = 42
        self.assertNotEqual(copy, self.mat1)
        self.assertEqual(copy[0,0], 42)

        # zero element -> nonzero element
        copy = self.mat1.copy()
        copy[0,1] = 42
        self.assertNotEqual(copy, self.mat1)
        self.assertEqual(copy[0,1], 42)

    def test_getitem(self):
        """Test getting an element from the matrix."""
        with self.assertRaises(IndexError):
            _ = self.null[0,0]

        with self.assertRaises(IndexError):
            _ = self.empty[0]

        with self.assertRaises(IndexError):
            _ = self.empty[:,:]

        with self.assertRaises(IndexError):
            _ = self.empty[0:1,0]

        with self.assertRaises(IndexError):
            _ = self.empty[0,0:1]

        exp = ScipySparseMat(2,1)
        obs = self.empty[:,0]
        self.assertEqual(obs, exp)

        # Extracting a column.
        obs = self.mat1[:,2]
        self.assertEqual(obs, self.mat1.getCol(2))

        # Extracting a row.
        obs = self.mat1[1,:]
        self.assertEqual(obs, self.mat1.getRow(1))

        # Extracting a single element.
        self.assertEqual(self.empty[1,1], 0)
        self.assertEqual(self.mat1[1,2], 4)

        with self.assertRaises(IndexError):
            _ = self.mat1[1,3]


if __name__ == '__main__':
    main()
