#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from numpy import array
from scipy.sparse import lil_matrix
from biom.unit_test import TestCase, main
from biom.backends.scipysparse import (ScipySparseMat, to_scipy,
                                       list_nparray_to_scipy,
                                       list_list_to_scipy, list_scipy_to_scipy,
                                       nparray_to_scipy, dict_to_scipy,
                                       list_dict_to_scipy, coo_arrays_to_scipy)

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"


class ScipySparseMatTests(TestCase):

    def setUp(self):
        # 1 0 2
        # 3 0 4
        self.mat1 = ScipySparseMat(2, 3, data=array([[1, 0, 2], [3, 0, 4]]))

        # Empty/null cases (i.e., 0x0, 0xn, nx0).
        self.null1 = ScipySparseMat(0, 0)
        self.null2 = ScipySparseMat(0, 42)
        self.null3 = ScipySparseMat(42, 0)
        self.nulls = [self.null1, self.null2, self.null3]

        # 0 0
        # 0 0
        self.empty = ScipySparseMat(2, 2)

        # 1 0 3
        self.row_vec = ScipySparseMat(1, 3)
        self.row_vec[0, 0] = 1
        self.row_vec[0, 2] = 3

        # 1
        # 0
        # 3
        self.col_vec = ScipySparseMat(3, 1)
        self.col_vec[0, 0] = 1
        self.col_vec[2, 0] = 3

        # 1x1
        self.single_ele = ScipySparseMat(1, 1)
        self.single_ele[0, 0] = 42

        # Explicit zeros.
        self.explicit_zeros = ScipySparseMat(2, 3,
                                             data=([1, 2, 3, 0, 4], ([0, 0, 1, 1, 1], [0, 2, 0, 1, 2])))

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
        """Differentiate empty matrix from non-empty matrix."""
        for m in self.nulls:
            self.assertTrue(m.is_empty)

        self.assertFalse(self.mat1.is_empty)

    def test_shape(self):
        """What kind of shape are you in?"""
        self.assertEqual(self.null1.shape, (0, 0))
        self.assertEqual(self.null2.shape, (0, 42))
        self.assertEqual(self.null3.shape, (42, 0))
        self.assertEqual(self.mat1.shape, (2, 3))
        self.assertEqual(self.empty.shape, (2, 2))
        self.assertEqual(self.row_vec.shape, (1, 3))
        self.assertEqual(self.col_vec.shape, (3, 1))
        self.assertEqual(self.single_ele.shape, (1, 1))

    def test_dtype(self):
        """What's your type?"""
        for m in self.nulls:
            self.assertEqual(m.dtype, None)

        self.assertEqual(self.empty.dtype, float)
        self.assertEqual(self.row_vec.dtype, float)

    def test_fmt(self):
        """What format are you in?"""
        for m in self.nulls:
            self.assertEqual(m.fmt, None)

        self.assertEqual(self.empty.fmt, 'coo')
        self.assertEqual(self.mat1.fmt, 'csr')
        self.assertEqual(self.single_ele.fmt, 'lil')

    def test_size(self):
        """What is your NNZ?"""
        for m in self.nulls:
            self.assertEqual(m.size, 0)

        self.assertEqual(self.empty.size, 0)
        self.assertEqual(self.single_ele.size, 1)
        self.assertEqual(self.mat1.size, 4)
        self.assertEqual(self.explicit_zeros.size, 4)

    def test_convert(self):
        """Test sparse format conversion."""
        self.assertEqual(self.mat1.fmt, 'csr')
        self.mat1.convert('coo')
        self.assertEqual(self.mat1.fmt, 'coo')

        for m in self.nulls:
            self.assertEqual(m.fmt, None)
            m.convert('coo')
            self.assertEqual(m.fmt, None)

    def test_transpose(self):
        """Test transposition."""
        obs = self.null1.T
        self.assertEqual(obs, self.null1)
        self.assertFalse(obs is self.null1)

        obs = self.null2.T
        self.assertEqual(obs, self.null3)
        self.assertFalse(obs is self.null3)

        obs = self.null3.T
        self.assertEqual(obs, self.null2)
        self.assertFalse(obs is self.null2)

        obs = self.single_ele.T
        self.assertEqual(obs, self.single_ele)
        self.assertFalse(obs is self.single_ele)

        exp = ScipySparseMat(3, 2, data=array([[1, 3], [0, 0], [2, 4]]))
        obs = self.mat1.T
        self.assertEqual(obs, exp)

    def test_sum(self):
        """Test summing a matrix."""
        for m in self.nulls:
            self.assertEqual(m.sum(), 0)

        self.assertEqual(self.mat1.sum(), 10)
        self.assertEqual(self.mat1.sum(0), array([4, 0, 6]))
        self.assertEqual(self.mat1.sum(1), array([3, 7]))
        self.assertEqual(self.row_vec.sum(1), array([4]))
        self.assertEqual(self.col_vec.sum(0), array([4]))
        with self.assertRaises(ValueError):
            _ = self.mat1.sum(3)

    def test_getRow(self):
        """Test grabbing a row from the matrix."""
        for m in self.nulls:
            with self.assertRaises(IndexError):
                _ = m.getRow(0)

        exp = ScipySparseMat(1, 3, data=array([[1, 0, 2]]))
        obs = self.mat1.getRow(0)
        self.assertEqual(obs, exp)

        self.assertEqual(self.row_vec.getRow(0), self.row_vec)

    def test_getCol(self):
        """Test grabbing a column from the matrix."""
        for m in self.nulls:
            with self.assertRaises(IndexError):
                _ = m.getCol(0)

        exp = ScipySparseMat(2, 1, data=array([[1], [3]]))
        obs = self.mat1.getCol(0)
        self.assertEqual(obs, exp)

        self.assertEqual(self.col_vec.getCol(0), self.col_vec)

    def test_items_iteritems(self):
        """Test getting a list of non-zero elements."""
        exp = []
        for m in self.nulls + [self.empty]:
            self.assertEqual(m.items(), exp)
            self.assertEqual(list(m.iteritems()), exp)

        exp = [((0, 0), 1), ((0, 2), 2), ((1, 0), 3), ((1, 2), 4)]
        self.assertEqual(sorted(self.mat1.items()), exp)
        self.assertEqual(sorted(self.mat1.iteritems()), exp)

    def test_copy(self):
        """Test copying the matrix."""
        for m in self.nulls:
            copy = m.copy()
            self.assertEqual(copy, m)
            self.assertFalse(copy is m)

        copy = self.mat1.copy()
        self.assertEqual(copy, self.mat1)
        self.assertFalse(copy is self.mat1)

        copy[1, 1] = 42
        self.assertNotEqual(copy, self.mat1)

    def test_eq(self):
        """Test whether two matrices are equal."""
        self.assertTrue(self.null1 == ScipySparseMat(0, 0))
        self.assertTrue(self.null2 == ScipySparseMat(0, 42))
        self.assertTrue(self.null3 == ScipySparseMat(42, 0))
        self.assertTrue(self.empty == ScipySparseMat(2, 2))

        mat2 = ScipySparseMat(2, 3, data=array([[1, 0, 2], [3, 0, 4]]))
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
        self.assertTrue(self.null1 != array([]))

        # Wrong shape.
        self.assertTrue(self.null2 != self.null3)
        self.assertTrue(self.empty != ScipySparseMat(2, 1))

        # Wrong dtype.
        self.assertTrue(self.empty != ScipySparseMat(2, 2, dtype=int))

        # Wrong size.
        wrong_size = ScipySparseMat(2, 2)
        self.assertTrue(self.empty == wrong_size)
        wrong_size[1, 0] = 42
        self.assertTrue(self.empty != wrong_size)

        # Wrong size.
        wrong_data = self.mat1.copy()
        self.assertTrue(self.mat1 == wrong_data)
        wrong_data[0, 2] = 42
        self.assertTrue(self.mat1 != wrong_data)
        self.assertTrue(wrong_data != self.mat1)

    def test_str(self):
        """Test getting string representation of the matrix."""
        for m in self.nulls:
            self.assertEqual(str(m), '<%dx%d empty/null sparse matrix>' %
                             (m.shape[0], m.shape[1]))

    def test_setitem(self):
        """Test setting an element in the matrix."""
        for m in self.nulls:
            with self.assertRaises(IndexError):
                m[0, 0] = 42

        with self.assertRaises(IndexError):
            self.empty[0] = [42, 42]

        with self.assertRaises(ValueError):
            self.mat1[0, 0] = 0

        with self.assertRaises(ValueError):
            self.mat1[0, 0] = 0.0

        # Setting existing zero element doesn't change matrix.
        copy = self.mat1.copy()
        copy[0, 1] = 0.0
        self.assertEqual(copy, self.mat1)

        # Setting existing nonzero element to the same thing doesn't change
        # matrix.
        copy = self.mat1.copy()
        copy[0, 0] = 1.0
        self.assertEqual(copy, self.mat1)

        # nonzero element -> nonzero element
        copy = self.mat1.copy()
        copy[0, 0] = 42
        self.assertNotEqual(copy, self.mat1)
        self.assertEqual(copy[0, 0], 42)

        # zero element -> nonzero element
        copy = self.mat1.copy()
        copy[0, 1] = 42
        self.assertNotEqual(copy, self.mat1)
        self.assertEqual(copy[0, 1], 42)

    def test_getitem(self):
        """Test getting an element from the matrix."""
        for m in self.nulls:
            with self.assertRaises(IndexError):
                _ = m[0, 0]

        with self.assertRaises(IndexError):
            _ = self.empty[0]

        with self.assertRaises(IndexError):
            _ = self.empty[:, :]

        with self.assertRaises(IndexError):
            _ = self.empty[0:1, 0]

        with self.assertRaises(IndexError):
            _ = self.empty[0, 0:1]

        exp = ScipySparseMat(2, 1)
        obs = self.empty[:, 0]
        self.assertEqual(obs, exp)

        # Extracting a column.
        obs = self.mat1[:, 2]
        self.assertEqual(obs, self.mat1.getCol(2))

        # Extracting a row.
        obs = self.mat1[1, :]
        self.assertEqual(obs, self.mat1.getRow(1))

        # Extracting a single element.
        self.assertEqual(self.empty[1, 1], 0)
        self.assertEqual(self.mat1[1, 2], 4)

        with self.assertRaises(IndexError):
            _ = self.mat1[1, 3]

# These tests are pretty much copied from CSMat's conversion tests...


class SupportTests(TestCase):

    def test_coo_arrays_to_scipy(self):
        """convert (values, (row, col)) to scipy"""
        n_rows, n_cols = 3, 4
        exp_d = lil_matrix((n_rows, n_cols))
        exp_d[0, 0] = 10
        exp_d[1, 3] = 5
        exp_d[2, 1] = 2
        exp_d = exp_d.tocoo()
        exp = ScipySparseMat(n_rows, n_cols, data=exp_d)

        data = (array([5.0, 2.0, 10.0]), (array([1, 2, 0]), array([3, 1, 0])))
        obs = coo_arrays_to_scipy(data, shape=(n_rows, n_cols))
        self.assertEqual(obs, exp)

    def test_list_list_to_scipy(self):
        """convert [[row,col,value], ...] to scipy"""
        input = [[0, 0, 1], [1, 1, 5.0], [0, 2, 6]]
        exp = ScipySparseMat(2, 3)
        exp[0, 0] = 1.0
        exp[1, 1] = 5.0
        exp[0, 2] = 6
        obs = list_list_to_scipy(input)
        self.assertEqual(obs, exp)

    def test_nparray_to_scipy(self):
        """Convert nparray to scipy"""
        input = array([[1, 2, 3, 4], [-1, 6, 7, 8], [9, 10, 11, 12]])
        exp = ScipySparseMat(3, 4)
        exp[0, 0] = 1
        exp[0, 1] = 2
        exp[0, 2] = 3
        exp[0, 3] = 4
        exp[1, 0] = -1
        exp[1, 1] = 6
        exp[1, 2] = 7
        exp[1, 3] = 8
        exp[2, 0] = 9
        exp[2, 1] = 10
        exp[2, 2] = 11
        exp[2, 3] = 12
        obs = nparray_to_scipy(input)
        self.assertEqual(obs, exp)

    def test_list_dict_to_scipy(self):
        """Take a list of dicts and condense down to a single dict"""
        input = [{(0, 0): 10, (0, 1): 2}, {(1, 2): 15}, {(0, 3): 7}]
        exp = ScipySparseMat(3, 4)
        exp[0, 0] = 10
        exp[0, 1] = 2
        exp[1, 2] = 15
        exp[2, 3] = 7
        obs = list_dict_to_scipy(input)
        self.assertEqual(obs, exp)

    def test_dict_to_scipy(self):
        """Take a dict and convert to scipy"""
        input = {(0, 1): 5, (1, 0): 2, (2, 1): 6}
        exp = ScipySparseMat(3, 2)
        exp[(0, 1)] = 5
        exp[(1, 0)] = 2
        exp[(2, 1)] = 6
        obs = dict_to_scipy(input)
        self.assertEqual(obs, exp)

    def test_to_scipy(self):
        """Convert to expected scipy types"""
        vals = {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8}
        obs = to_scipy(vals)
        exp = ScipySparseMat(2, 2)
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(1, 0)] = 7
        exp[(1, 1)] = 8
        self.assertEqual(obs, exp)

        input = {(0, 1): 5, (10, 8): -1.23}

        exp = ScipySparseMat(11, 9)
        exp[(0, 1)] = 5
        exp[(10, 8)] = -1.23
        obs = to_scipy(input)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # test transpose
        exp = ScipySparseMat(9, 11)
        exp[(1, 0)] = 5
        exp[(8, 10)] = -1.23
        obs = to_scipy(input, transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of dicts, transpose
        exp = ScipySparseMat(3, 2)
        exp[(0, 0)] = 5.0
        exp[(1, 0)] = 6.0
        exp[(2, 0)] = 7.0
        exp[(0, 1)] = 8.0
        exp[(1, 1)] = 9.0
        exp[(2, 1)] = 10.0
        obs = to_scipy([{(0, 0): 5, (0, 1): 6, (0, 2): 7},
                        {(1, 0): 8, (1, 1): 9, (1, 2): 10}],
                       transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of ScipySparseMats
        exp = ScipySparseMat(2, 3)
        exp[(0, 0)] = 5
        exp[(0, 1)] = 6
        exp[(0, 2)] = 7
        exp[(1, 0)] = 8
        exp[(1, 1)] = 9
        exp[(1, 2)] = 10
        row1 = ScipySparseMat(1, 3)
        row1[(0, 0)] = 5
        row1[(0, 1)] = 6
        row1[(0, 2)] = 7
        row2 = ScipySparseMat(1, 3)
        row2[(0, 0)] = 8
        row2[(0, 1)] = 9
        row2[(0, 2)] = 10
        obs = to_scipy([row1, row2])
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # test empty set
        exp = ScipySparseMat(0, 0)
        obs = to_scipy([])
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

    def test_list_nparray_to_scipy(self):
        """lists of nparrays to scipy"""
        ins = [array([0, 2, 1, 0]), array([1, 0, 0, 1])]
        exp = ScipySparseMat(2, 4)
        exp[0, 1] = 2
        exp[0, 2] = 1
        exp[1, 0] = 1
        exp[1, 3] = 1
        obs = list_nparray_to_scipy(ins)
        self.assertEqual(obs, exp)

    def test_list_scipy_to_scipy(self):
        """list of ScipySparseMats to ScipySparseMat"""
        ins = [ScipySparseMat(1, 4), ScipySparseMat(1, 4)]
        ins[0][0, 0] = 5
        ins[0][0, 1] = 10
        ins[1][0, 2] = 1
        ins[1][0, 3] = 2
        exp = ScipySparseMat(2, 4)
        exp[0, 0] = 5
        exp[0, 1] = 10
        exp[1, 2] = 1
        exp[1, 3] = 2
        obs = list_scipy_to_scipy(ins)
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
