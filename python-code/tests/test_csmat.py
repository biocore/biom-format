#!/usr/bin/env python

from numpy import zeros, ndarray, array
from biom.unit_test import TestCase, main
from biom.table import flatten
from biom.csmat import CSMat, to_csmat, \
    list_nparray_to_csmat, list_list_to_csmat, \
    list_csmat_to_csmat, nparray_to_csmat, \
    dict_to_csmat, list_dict_to_csmat

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
   
class CSMatTests(TestCase):
    def setUp(self):
        # No empty rows/columns.

        # 1 2 0 0
        # 0 0 3 0
        # 0 0 0 4
        self.obj = CSMat(3,4)
        self.obj._coo_rows = [0, 0, 1, 2]
        self.obj._coo_cols = [0, 1, 2, 3]
        self.obj._coo_values = [1, 2, 3, 4]

        # Empty rows.

        # 0 0 0 0
        # 1 9 0 0
        # 1 1 2 0
        self.empty_row_start = CSMat(3, 4)
        self.empty_row_start._coo_values = [1, 9, 1, 1, 2]
        self.empty_row_start._coo_rows = [1, 1, 2, 2, 2]
        self.empty_row_start._coo_cols = [0, 1, 0, 1, 2]

        # 1 9 0 0
        # 0 0 0 0
        # 1 1 2 0
        self.empty_row_mid = CSMat(3, 4)
        self.empty_row_mid._coo_values = [1, 9, 1, 1, 2]
        self.empty_row_mid._coo_rows = [0, 0, 2, 2, 2]
        self.empty_row_mid._coo_cols = [0, 1, 0, 1, 2]

        # 1 9 0 0
        # 1 1 2 0
        # 0 0 0 0
        self.empty_row_end = CSMat(3, 4)
        self.empty_row_end._coo_values = [1, 9, 1, 1, 2]
        self.empty_row_end._coo_rows = [0, 0, 1, 1, 1]
        self.empty_row_end._coo_cols = [0, 1, 0, 1, 2]

        # 1 9 0 0
        # 0 0 0 0
        # 0 0 0 0
        self.empty_rows = CSMat(3, 4)
        self.empty_rows._coo_values = [1, 9]
        self.empty_rows._coo_rows = [0, 0]
        self.empty_rows._coo_cols = [0, 1]

        # Empty columns.

        # 0 0 0
        # 0 9 1
        # 0 8 0
        # 0 0 3
        self.empty_col_start = CSMat(4, 3)
        self.empty_col_start._coo_values = [9, 1, 8, 3]
        self.empty_col_start._coo_rows = [1, 1, 2, 3]
        self.empty_col_start._coo_cols = [1, 2, 1, 2]

        # 0 0 0
        # 9 0 1
        # 8 0 0
        # 0 0 3
        self.empty_col_mid = CSMat(4, 3)
        self.empty_col_mid._coo_values = [9, 1, 8, 3]
        self.empty_col_mid._coo_rows = [1, 1, 2, 3]
        self.empty_col_mid._coo_cols = [0, 2, 0, 2]

        # 0 0 0
        # 9 1 0
        # 8 0 0
        # 0 3 0
        self.empty_col_end = CSMat(4, 3)
        self.empty_col_end._coo_values = [9, 1, 8, 3]
        self.empty_col_end._coo_rows = [1, 1, 2, 3]
        self.empty_col_end._coo_cols = [0, 1, 0, 1]

        # 0 0 0
        # 0 0 1
        # 0 0 0
        # 0 0 3
        self.empty_cols = CSMat(4, 3)
        self.empty_cols._coo_values = [1, 3]
        self.empty_cols._coo_rows = [1, 3]
        self.empty_cols._coo_cols = [2, 2]

        # Empty matrix.

        # 0 0 0 0
        # 0 0 0 0
        # 0 0 0 0
        self.empty = CSMat(3, 4)

    def test_copy(self):
        """copy thy self"""
        obs = self.obj.copy()
        self.assertEqual(obs, self.obj)
        self.assertNotEqual(id(obs), id(self.obj))
        self.assertNotEqual(id(obs._coo_rows), id(self.obj._coo_rows))
        self.assertNotEqual(id(obs._coo_cols), id(self.obj._coo_cols))
        self.assertNotEqual(id(obs._coo_values), id(self.obj._coo_values))
        self.assertNotEqual(id(obs._pkd_ax), id(self.obj._pkd_ax))
        self.assertNotEqual(id(obs._unpkd_ax), id(self.obj._unpkd_ax))
        self.assertNotEqual(id(obs._values), id(self.obj._values))
        obs[1,2] = 10
        self.assertNotEqual(obs, self.obj)
        obs[1,2] = 3
        obs[0,0] = 5
        self.assertNotEqual(obs, self.obj)

    def test_absorbUpdates(self):
        """absorb all the updates"""
        self.obj.convert("csr")
        self.obj.bulkCOOUpdate([0,2],[3,0],[-10,50])
        self.obj.absorbUpdates()
        self.assertEqual(self.obj._order, "csr")
        self.assertEqual(self.obj._pkd_ax, [0,3,4,6])
        self.assertEqual(self.obj._unpkd_ax, [0,1,3,2,0,3])
        self.assertEqual(self.obj._values, [1,2,-10,3,50,4])

    def test_setitem(self):
        self.obj[(2,2)] = 10
        exp = sorted([((0,0),1),((0,1),2),((2,2),10),((1,2),3),((2,3),4)])
        self.assertEqual(sorted(self.obj.items()), exp)
        self.assertRaises(IndexError, self.obj.__setitem__, (100,50), 10)
        self.assertRaises(ValueError, self.empty_cols.__setitem__, [3,2], 0)

        self.empty_cols.convert("csr")
        self.empty_cols[(2,2)] = 42
        exp = sorted([((1,2),1),((2,2),42),((3,2),3)])
        self.assertEqual(sorted(self.empty_cols.items()), exp)

    def test_getitem_simple(self):
        """Tests simple getitem"""
        self.assertEqual(self.obj[(1,2)], 3)
        self.assertEqual(self.obj[1,2], 3)
        self.assertEqual(self.obj[1,1], 0)
        self.assertRaises(IndexError, self.obj.__getitem__, (3,3))
        self.assertRaises(IndexError, self.obj.__getitem__, (-1,2))
        self.assertRaises(IndexError, self.obj.__getitem__, 1)

    def test_getitem_slice(self):
        """Tests for slices on getitem"""
        exp = CSMat(1,4)
        exp[0,0] = 1
        exp[0,1] = 2
        self.assertEqual(exp, self.obj[0,:])

        exp = CSMat(3,1)
        exp[1,0] = 3
        self.assertEqual(exp, self.obj[:,2])

        self.assertRaises(IndexError, self.obj.__getitem__, (10,slice(None)))
        self.assertRaises(AttributeError, self.obj.__getitem__, (3, slice(1,2,3)))

        exp = CSMat(4,1)
        self.assertEqual(self.empty_cols[:,0], exp)
        self.assertEqual(self.empty_cols[:,1], exp)
        exp[1,0] = 1
        exp[3,0] = 3
        self.assertEqual(self.empty_cols[:,2], exp)

        exp = CSMat(1,3)
        self.assertEqual(self.empty_cols[0,:], exp)
        self.assertEqual(self.empty_cols[2,:], exp)
        exp[0,2] = 1
        self.assertEqual(self.empty_cols[1,:], exp)
        exp = CSMat(1,3)
        exp[0,2] = 3
        self.assertEqual(self.empty_cols[3,:], exp)

    def test_getitem_direct(self):
        """test _getitem"""
        self.assertEqual(self.obj._order, "coo")
        self.assertEqual(self.obj._getitem((0,0)), (0, 0, 0))
        self.assertEqual(self.obj._getitem((0,1)), (1, 1, 1))
        self.assertEqual(self.obj._getitem((1,2)), (2, 2, 2))
        self.assertEqual(self.obj._getitem((2,3)), (3, 3, 3))
        self.assertEqual(self.obj._getitem((1,1)), (None, None, None))
        self.assertEqual(self.obj._order, "coo")

        self.obj.convert("csr")
        self.assertEqual(self.obj._order, "csr")
        self.assertEqual(self.obj._getitem((0,0)), (0, 0, 0))
        self.assertEqual(self.obj._getitem((0,1)), (0, 1, 1))
        self.assertEqual(self.obj._getitem((1,2)), (1, 2, 2))
        self.assertEqual(self.obj._getitem((2,3)), (2, 3, 3))
        self.assertEqual(self.obj._getitem((1,1)), (None, None, None))
        self.assertEqual(self.obj._order, "csr")

        self.obj.convert("csc")
        self.assertEqual(self.obj._order, "csc")
        self.assertEqual(self.obj._getitem((0,0)), (0, 0, 0))
        self.assertEqual(self.obj._getitem((0,1)), (1, 1, 1))
        self.assertEqual(self.obj._getitem((1,2)), (2, 2, 2))
        self.assertEqual(self.obj._getitem((2,3)), (3, 3, 3))
        self.assertEqual(self.obj._getitem((1,1)), (None, None, None))
        self.assertEqual(self.obj._order, "csc")

        self.assertEqual(self.empty_cols._order, "coo")
        self.assertEqual(self.empty_cols._getitem((2,2)), (None, None, None))
        self.assertEqual(self.empty_cols._getitem((1,2)), (0, 0, 0))

        self.empty_cols.convert("csr")
        self.assertEqual(self.empty_cols._order, "csr")
        self.assertEqual(self.empty_cols._getitem((2,2)), (None, None, None))
        self.assertEqual(self.empty_cols._getitem((1,2)), (1, 0, 0))

        self.empty_cols.convert("csc")
        self.assertEqual(self.empty_cols._order, "csc")
        self.assertEqual(self.empty_cols._getitem((2,2)), (None, None, None))
        self.assertEqual(self.empty_cols._getitem((1,2)), (0, 2, 0))

        self.assertEqual(self.empty._order, "coo")
        self.assertEqual(self.empty._getitem((0,0)), (None, None, None))
        self.assertEqual(self.empty._getitem((2,2)), (None, None, None))

        self.empty.convert("csr")
        self.assertEqual(self.empty._getitem((0,0)), (None, None, None))
        self.assertEqual(self.empty._getitem((2,2)), (None, None, None))

        self.empty.convert("csc")
        self.assertEqual(self.empty._getitem((0,0)), (None, None, None))
        self.assertEqual(self.empty._getitem((2,2)), (None, None, None))

    def test_items(self):
        """Get items out"""
        exp = sorted([((0,0),1.0),((0,1),2.0),((1,2),3.0),((2,3),4.0)])
        obs_coo = sorted(self.obj.items())
        self.obj.convert("csr")
        obs_csr = sorted(self.obj.items())
        self.obj.convert("csc")
        obs_csc = sorted(self.obj.items())

        self.assertEqual(obs_coo, exp)
        self.assertEqual(obs_csr, exp)
        self.assertEqual(obs_csc, exp)

        exp = sorted([((1,2),1),((3,2),3)])
        obs_coo = sorted(self.empty_cols.items())
        self.empty_cols.convert("csr")
        obs_csr = sorted(self.empty_cols.items())
        self.empty_cols.convert("csc")
        obs_csc = sorted(self.empty_cols.items())

        self.assertEqual(obs_coo, exp)
        self.assertEqual(obs_csr, exp)
        self.assertEqual(obs_csc, exp)

        exp = []
        obs_coo = sorted(self.empty.items())
        self.empty.convert("csr")
        obs_csr = sorted(self.empty.items())
        self.empty.convert("csc")
        obs_csc = sorted(self.empty.items())

        self.assertEqual(obs_coo, exp)
        self.assertEqual(obs_csr, exp)
        self.assertEqual(obs_csc, exp)

    def test_iteritems(self):
        """Get items out"""
        exp = sorted([((0,0),1.0),((0,1),2.0),((1,2),3.0),((2,3),4.0)])
        obs_coo = sorted(self.obj.iteritems())
        self.obj.convert("csr")
        obs_csr = sorted(self.obj.iteritems())
        self.obj.convert("csc")
        obs_csc = sorted(self.obj.iteritems())

        self.assertEqual(obs_coo, exp)
        self.assertEqual(obs_csr, exp)
        self.assertEqual(obs_csc, exp)

        exp = sorted([((1,2),1),((3,2),3)])
        obs_coo = sorted(self.empty_cols.iteritems())
        self.empty_cols.convert("csr")
        obs_csr = sorted(self.empty_cols.iteritems())
        self.empty_cols.convert("csc")
        obs_csc = sorted(self.empty_cols.iteritems())

        self.assertEqual(obs_coo, exp)
        self.assertEqual(obs_csr, exp)
        self.assertEqual(obs_csc, exp)

        exp = []
        obs_coo = sorted(self.empty.iteritems())
        self.empty.convert("csr")
        obs_csr = sorted(self.empty.iteritems())
        self.empty.convert("csc")
        obs_csc = sorted(self.empty.iteritems())

        self.assertEqual(obs_coo, exp)
        self.assertEqual(obs_csr, exp)
        self.assertEqual(obs_csc, exp)

    def test_contains(self):
        """Make sure we can check things exist"""
        sm1 = CSMat(3,4)
        for r in range(3):
            for c in range(4):
                assert (r,c) not in sm1
        sm1[1,2] = 0
        assert (1,2) not in sm1
        sm1[1,2] = 10
        assert (1,2) in sm1

        assert (0,2) not in self.empty_cols
        assert (1,2) in self.empty_cols
        self.empty_cols.convert("csr")
        assert (0,2) not in self.empty_cols
        assert (1,2) in self.empty_cols

    def test_eq(self):
        """Tests for equality"""
        sm1 = CSMat(2,3)
        sm1.bulkCOOUpdate([0,1,1],[0,1,2],[1,2,3])
        sm2 = CSMat(2,3)
        sm2.bulkCOOUpdate([0,1,1],[0,1,2],[1,2,3])
        sm3 = CSMat(3,2)
        sm3.bulkCOOUpdate([0,1,2],[0,1,1],[1,2,3])

        sm1.convert('csr')
        self.assertEqual(sm1, sm2)
        self.assertNotEqual(sm1,sm3)

        sm1[0,1] = 10
        sm2[0,1] = 5
        self.assertNotEqual(sm1, sm2)
        sm1.absorbUpdates()
        sm2.absorbUpdates()
        sm2[0,1] = 10
        self.assertEqual(sm1, sm2)

        # Test empty rows/columns.
        self.empty_cols.convert("csr")
        self.empty_cols.update({(1,0):2})
        exp = CSMat(4, 3)
        exp[1,0] = 2
        exp[1,2] = 1
        exp[3,2] = 3
        self.assertEqual(self.empty_cols, exp)

        # Test empty matrices.
        empty_sm1 = CSMat(2,3)
        empty_sm2 = CSMat(2,3)
        empty_sm3 = CSMat(2,2)
        self.assertTrue(empty_sm1 == empty_sm2)
        self.assertFalse(empty_sm1 == empty_sm3)

        empty_sm1.convert("csr")
        empty_sm2.convert("csc")
        self.assertTrue(empty_sm1 == empty_sm2)
        empty_sm2.convert("coo")
        self.assertTrue(empty_sm1 == empty_sm2)

    def test_getRow(self):
        """Get a row"""
        exp = CSMat(1,4)
        exp.update({(0,2):3})
        obs = self.obj.getRow(1)
        self.assertEqual(obs, exp)

        exp = CSMat(1,4)
        obs = self.obj.getRow(2)
        exp.update({(0,3):4})
        self.assertEqual(obs,exp)

        obj = CSMat(4,1)
        obj[0,0] = 5
        obj[1,0] = 6
        obj[2,0] = 7
        obj[3,0] = 8
        exp1 = CSMat(1,1)
        exp1[0,0] = 5
        exp2 = CSMat(1,1)
        exp2[0,0] = 6
        exp3 = CSMat(1,1)
        exp3[0,0] = 7
        exp4 = CSMat(1,1)
        exp4[0,0] = 8
        obs1 = obj.getRow(0)
        obs2 = obj.getRow(1)
        obs3 = obj.getRow(2)
        obs4 = obj.getRow(3)
        self.assertRaises(IndexError,obj.getRow, 4)

        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)
        self.assertEqual(obs4, exp4)

        self.assertRaises(IndexError, self.obj.getRow, -1)

        # Test matrices with empty rows.
        exp = CSMat(1, 4)
        exp[0,0] = 1
        exp[0,1] = 9
        obs = self.empty_row_start.getRow(1)
        self.assertEqual(obs, exp)

        exp = CSMat(1, 4)
        obs = self.empty_row_start.getRow(0)
        self.assertEqual(obs, exp)

        obs = self.empty_row_mid.getRow(1)
        self.assertEqual(obs, exp)

        obs = self.empty_row_end.getRow(2)
        self.assertEqual(obs, exp)

        obs = self.empty_rows.getRow(1)
        self.assertEqual(obs, exp)

        obs = self.empty_rows.getRow(2)
        self.assertEqual(obs, exp)

        # Test completely empty matrix.
        obs = self.empty.getRow(2)
        self.assertEqual(obs, exp)

    def test_getCol(self):
        """Get a col"""
        exp = CSMat(3,1)
        exp.update({(1,0):3.0})
        obs = self.obj.getCol(2)
        self.assertEqual(obs,exp)

        exp = CSMat(3,1)
        exp.update({(0,0):2})
        obs = self.obj.getCol(1)
        self.assertEqual(obs,exp)

        self.assertRaises(IndexError, self.obj.getCol, -1)

        # Test matrices with empty columns.
        exp = CSMat(4, 1)
        exp[1,0] = 9
        exp[2,0] = 8
        obs = self.empty_col_start.getCol(1)
        self.assertEqual(obs, exp)

        exp = CSMat(4, 1)
        obs = self.empty_col_start.getCol(0)
        self.assertEqual(obs, exp)

        obs = self.empty_col_mid.getCol(1)
        self.assertEqual(obs, exp)

        obs = self.empty_col_end.getCol(2)
        self.assertEqual(obs, exp)

        obs = self.empty_cols.getCol(0)
        self.assertEqual(obs, exp)

        obs = self.empty_cols.getCol(1)
        self.assertEqual(obs, exp)

        # Test completely empty matrix.
        exp = CSMat(3, 1)
        obs = self.empty.getCol(2)
        self.assertEqual(obs, exp)

    def test_update(self):
        """updates should work on new and inplace values"""
        self.obj.update({(0,0):10,(0,3):6})
        exp = CSMat(3,4)
        exp[0,0] = 10
        exp[0,1] = 2
        exp[0,3] = 6
        exp[1,2] = 3
        exp[2,3] = 4
        self.assertEqual(exp, self.obj)

        self.empty_cols.convert("csr")
        self.empty_cols.update({(1,0):2})
        exp = CSMat(4, 3)
        exp[1,0] = 2
        exp[1,2] = 1
        exp[3,2] = 3
        self.assertEqual(self.empty_cols, exp)

    def test_transpose(self):
        """test transpose"""
        exp = CSMat(4,3)
        exp.update({(0,0):1,(1,0):2,(2,1):3,(3,2):4})
        obs = self.obj.T
        self.assertEqual(obs, exp)

        exp = CSMat(3, 4)
        exp.update({(2,1):1,(2,3):3})
        obs = self.empty_cols.T
        self.assertEqual(obs, exp)

        exp.convert("csc")
        self.empty_cols.convert("csc")
        obs = self.empty_cols.T
        self.assertEqual(obs, exp)

        exp = CSMat(4, 3)
        obs = self.empty.T
        self.assertEqual(obs, exp)

        exp.convert("csr")
        self.empty.convert("csr")
        obs = self.empty.T
        self.assertEqual(obs, exp)

    def test_bulkCOOUpdate(self):
        """Stages data"""
        self.obj.convert("csr")
        self.obj.bulkCOOUpdate([1,2],[3,4],[5,6])
        self.assertEqual(self.obj._coo_values, [5,6])
        self.assertEqual(self.obj._coo_cols, [3,4])
        self.assertEqual(self.obj._coo_rows, [1,2])

        # Make sure zeros are ignored.
        self.empty_row_start.convert("csc")
        self.empty_row_start.bulkCOOUpdate([0,2],[2,3],[42.0,0.00])
        self.assertEqual(self.empty_row_start._coo_values, [42.0])
        self.assertEqual(self.empty_row_start._coo_cols, [2])
        self.assertEqual(self.empty_row_start._coo_rows, [0])

    def test_hasUpdates(self):
        """Do we have updates?"""
        self.obj.convert("csr")
        self.assertFalse(self.obj.hasUpdates())
        self.obj.bulkCOOUpdate([1,2],[3,4],[5,6])
        self.assertTrue(self.obj.hasUpdates())

    def test_get_size(self):
        """Number of nonzeros"""
        self.assertEqual(self.obj.size, 4)
        self.obj.bulkCOOUpdate([1,2],[3,4],[5,6])
        self.assertEqual(self.obj.size, 6)

        self.assertEqual(self.empty_rows.size, 2)
        self.empty_rows.convert("csr")
        self.assertEqual(self.empty_rows.size, 2)
        self.empty_rows.convert("csc")
        self.assertEqual(self.empty_rows.size, 2)

        self.assertEqual(self.empty.size, 0)
        self.empty_rows.convert("csr")
        self.assertEqual(self.empty.size, 0)
        self.empty_rows.convert("csc")
        self.assertEqual(self.empty.size, 0)

    def test_convert_coo_csr(self):
        """convert coo to csr"""
        self.obj.convert("csr")
        self.assertEqual(self.obj._order, "csr")
        self.assertEqual(self.obj._coo_values, [])
        self.assertEqual(self.obj._coo_rows, [])
        self.assertEqual(self.obj._coo_cols, [])

        self.assertEqual(self.obj._pkd_ax, [0, 2, 3, 4])
        self.assertEqual(self.obj._unpkd_ax, [0, 1, 2, 3])
        self.assertEqual(self.obj._values, [1,2,3,4])

        # Test empty rows.
        self.empty_row_start.convert("csr")

        self.assertEqual(self.empty_row_start._order, "csr")
        self.assertEqual(self.empty_row_start._coo_values, [])
        self.assertEqual(self.empty_row_start._coo_rows, [])
        self.assertEqual(self.empty_row_start._coo_cols, [])

        self.assertEqual(self.empty_row_start._pkd_ax, [0, 0, 2, 5])
        self.assertEqual(self.empty_row_start._unpkd_ax, [0, 1, 0, 1, 2])
        self.assertEqual(self.empty_row_start._values, [1, 9, 1, 1, 2])

        self.empty_row_mid.convert("csr")

        self.assertEqual(self.empty_row_mid._order, "csr")
        self.assertEqual(self.empty_row_mid._coo_values, [])
        self.assertEqual(self.empty_row_mid._coo_rows, [])
        self.assertEqual(self.empty_row_mid._coo_cols, [])

        self.assertEqual(self.empty_row_mid._pkd_ax, [0, 2, 2, 5])
        self.assertEqual(self.empty_row_mid._unpkd_ax, [0, 1, 0, 1, 2])
        self.assertEqual(self.empty_row_mid._values, [1, 9, 1, 1, 2])

        self.empty_row_end.convert("csr")

        self.assertEqual(self.empty_row_end._order, "csr")
        self.assertEqual(self.empty_row_end._coo_values, [])
        self.assertEqual(self.empty_row_end._coo_rows, [])
        self.assertEqual(self.empty_row_end._coo_cols, [])

        self.assertEqual(self.empty_row_end._pkd_ax, [0, 2, 5, 5])
        self.assertEqual(self.empty_row_end._unpkd_ax, [0, 1, 0, 1, 2])
        self.assertEqual(self.empty_row_end._values, [1, 9, 1, 1, 2])

        self.empty_rows.convert("csr")

        self.assertEqual(self.empty_rows._order, "csr")
        self.assertEqual(self.empty_rows._coo_values, [])
        self.assertEqual(self.empty_rows._coo_rows, [])
        self.assertEqual(self.empty_rows._coo_cols, [])

        self.assertEqual(self.empty_rows._pkd_ax, [0, 2, 2, 2])
        self.assertEqual(self.empty_rows._unpkd_ax, [0, 1])
        self.assertEqual(self.empty_rows._values, [1, 9])

        # Test an empty matrix.
        self.empty.convert("csr")

        self.assertEqual(self.empty._order, "csr")
        self.assertEqual(self.empty._coo_values, [])
        self.assertEqual(self.empty._coo_rows, [])
        self.assertEqual(self.empty._coo_cols, [])

        self.assertEqual(self.empty._pkd_ax, [0, 0, 0, 0])
        self.assertEqual(self.empty._unpkd_ax, [])
        self.assertEqual(self.empty._values, [])

    def test_convert_coo_csc(self):
        """convert coo to csc"""
        self.obj.convert("csc")
        self.assertEqual(self.obj._order, "csc")
        self.assertEqual(self.obj._coo_values, [])
        self.assertEqual(self.obj._coo_rows, [])
        self.assertEqual(self.obj._coo_cols, [])

        self.assertEqual(self.obj._pkd_ax, [0, 1, 2, 3, 4])
        self.assertEqual(self.obj._unpkd_ax, [0, 0, 1, 2])
        self.assertEqual(self.obj._values, [1,2,3,4])

        # Test empty columns.
        self.empty_col_start.convert("csc")

        self.assertEqual(self.empty_col_start._order, "csc")
        self.assertEqual(self.empty_col_start._coo_values, [])
        self.assertEqual(self.empty_col_start._coo_rows, [])
        self.assertEqual(self.empty_col_start._coo_cols, [])

        self.assertEqual(self.empty_col_start._pkd_ax, [0, 0, 2, 4])
        self.assertEqual(self.empty_col_start._unpkd_ax, [1, 2, 1, 3])
        self.assertEqual(self.empty_col_start._values, [9, 8, 1, 3])

        self.empty_col_mid.convert("csc")

        self.assertEqual(self.empty_col_mid._order, "csc")
        self.assertEqual(self.empty_col_mid._coo_values, [])
        self.assertEqual(self.empty_col_mid._coo_rows, [])
        self.assertEqual(self.empty_col_mid._coo_cols, [])

        self.assertEqual(self.empty_col_mid._pkd_ax, [0, 2, 2, 4])
        self.assertEqual(self.empty_col_mid._unpkd_ax, [1, 2, 1, 3])
        self.assertEqual(self.empty_col_mid._values, [9, 8, 1, 3])

        self.empty_col_end.convert("csc")

        self.assertEqual(self.empty_col_end._order, "csc")
        self.assertEqual(self.empty_col_end._coo_values, [])
        self.assertEqual(self.empty_col_end._coo_rows, [])
        self.assertEqual(self.empty_col_end._coo_cols, [])

        self.assertEqual(self.empty_col_end._pkd_ax, [0, 2, 4, 4])
        self.assertEqual(self.empty_col_end._unpkd_ax, [1, 2, 1, 3])
        self.assertEqual(self.empty_col_end._values, [9, 8, 1, 3])

        self.empty_cols.convert("csc")

        self.assertEqual(self.empty_cols._order, "csc")
        self.assertEqual(self.empty_cols._coo_values, [])
        self.assertEqual(self.empty_cols._coo_rows, [])
        self.assertEqual(self.empty_cols._coo_cols, [])

        self.assertEqual(self.empty_cols._pkd_ax, [0, 0, 0, 2])
        self.assertEqual(self.empty_cols._unpkd_ax, [1, 3])
        self.assertEqual(self.empty_cols._values, [1, 3])

        # Test an empty matrix.
        self.empty.convert("csc")

        self.assertEqual(self.empty._order, "csc")
        self.assertEqual(self.empty._coo_values, [])
        self.assertEqual(self.empty._coo_rows, [])
        self.assertEqual(self.empty._coo_cols, [])

        self.assertEqual(self.empty._pkd_ax, [0, 0, 0, 0, 0])
        self.assertEqual(self.empty._unpkd_ax, [])
        self.assertEqual(self.empty._values, [])

    def test_convert_csr_coo(self):
        """convert csr to coo"""
        self.obj.convert("csr")
        self.obj.convert("coo")
        self.assertEqual(self.obj._order, "coo")
        self.assertEqual(self.obj._coo_values, [1,2,3,4])
        self.assertEqual(self.obj._coo_rows, [0,0,1,2])
        self.assertEqual(self.obj._coo_cols, [0,1,2,3])

        self.assertEqual(self.obj._pkd_ax, array([]))
        self.assertEqual(self.obj._unpkd_ax, array([]))
        self.assertEqual(self.obj._values, array([]))

        # Test empty rows.
        self.empty_row_start.convert("csr")
        self.empty_row_start.convert("coo")

        self.assertEqual(self.empty_row_start._order, "coo")
        self.assertEqual(self.empty_row_start._coo_values, [1,9,1,1,2])
        self.assertEqual(self.empty_row_start._coo_rows, [1,1,2,2,2])
        self.assertEqual(self.empty_row_start._coo_cols, [0,1,0,1,2])

        self.assertEqual(self.empty_row_start._pkd_ax, array([]))
        self.assertEqual(self.empty_row_start._unpkd_ax, array([]))
        self.assertEqual(self.empty_row_start._values, array([]))

        self.empty_row_mid.convert("csr")
        self.empty_row_mid.convert("coo")

        self.assertEqual(self.empty_row_mid._order, "coo")
        self.assertEqual(self.empty_row_mid._coo_values, [1,9,1,1,2])
        self.assertEqual(self.empty_row_mid._coo_rows, [0,0,2,2,2])
        self.assertEqual(self.empty_row_mid._coo_cols, [0,1,0,1,2])

        self.assertEqual(self.empty_row_mid._pkd_ax, array([]))
        self.assertEqual(self.empty_row_mid._unpkd_ax, array([]))
        self.assertEqual(self.empty_row_mid._values, array([]))

        self.empty_row_end.convert("csr")
        self.empty_row_end.convert("coo")

        self.assertEqual(self.empty_row_end._order, "coo")
        self.assertEqual(self.empty_row_end._coo_values, [1,9,1,1,2])
        self.assertEqual(self.empty_row_end._coo_rows, [0,0,1,1,1])
        self.assertEqual(self.empty_row_end._coo_cols, [0,1,0,1,2])

        self.assertEqual(self.empty_row_end._pkd_ax, array([]))
        self.assertEqual(self.empty_row_end._unpkd_ax, array([]))
        self.assertEqual(self.empty_row_end._values, array([]))

        self.empty_rows.convert("csr")
        self.empty_rows.convert("coo")

        self.assertEqual(self.empty_rows._order, "coo")
        self.assertEqual(self.empty_rows._coo_values, [1,9])
        self.assertEqual(self.empty_rows._coo_rows, [0,0])
        self.assertEqual(self.empty_rows._coo_cols, [0,1])

        self.assertEqual(self.empty_rows._pkd_ax, array([]))
        self.assertEqual(self.empty_rows._unpkd_ax, array([]))
        self.assertEqual(self.empty_rows._values, array([]))

        # Test an empty matrix.
        self.empty.convert("csr")
        self.empty.convert("coo")

        self.assertEqual(self.empty._order, "coo")
        self.assertEqual(self.empty._coo_values, [])
        self.assertEqual(self.empty._coo_rows, [])
        self.assertEqual(self.empty._coo_cols, [])

        self.assertEqual(self.empty._pkd_ax, array([]))
        self.assertEqual(self.empty._unpkd_ax, array([]))
        self.assertEqual(self.empty._values, array([]))

    def test_convert_csc_coo(self):
        """convert csc to coo"""
        self.obj.convert("csc")
        self.obj.convert("coo")
        self.assertEqual(self.obj._order, "coo")
        self.assertEqual(self.obj._coo_values, [1,2,3,4])
        self.assertEqual(self.obj._coo_rows, [0,0,1,2])
        self.assertEqual(self.obj._coo_cols, [0,1,2,3])

        self.assertEqual(self.obj._pkd_ax, array([]))
        self.assertEqual(self.obj._unpkd_ax, array([]))
        self.assertEqual(self.obj._values, array([]))

        # Test empty columns.
        self.empty_col_start.convert("csc")
        self.empty_col_start.convert("coo")

        self.assertEqual(self.empty_col_start._order, "coo")
        self.assertEqual(self.empty_col_start._coo_values, [9,8,1,3])
        self.assertEqual(self.empty_col_start._coo_rows, [1,2,1,3])
        self.assertEqual(self.empty_col_start._coo_cols, [1,1,2,2])

        self.assertEqual(self.empty_col_start._pkd_ax, array([]))
        self.assertEqual(self.empty_col_start._unpkd_ax, array([]))
        self.assertEqual(self.empty_col_start._values, array([]))

        self.empty_col_mid.convert("csc")
        self.empty_col_mid.convert("coo")

        self.assertEqual(self.empty_col_mid._order, "coo")
        self.assertEqual(self.empty_col_mid._coo_values, [9,8,1,3])
        self.assertEqual(self.empty_col_mid._coo_rows, [1,2,1,3])
        self.assertEqual(self.empty_col_mid._coo_cols, [0,0,2,2])

        self.assertEqual(self.empty_col_mid._pkd_ax, array([]))
        self.assertEqual(self.empty_col_mid._unpkd_ax, array([]))
        self.assertEqual(self.empty_col_mid._values, array([]))

        self.empty_col_end.convert("csc")
        self.empty_col_end.convert("coo")

        self.assertEqual(self.empty_col_end._order, "coo")
        self.assertEqual(self.empty_col_end._coo_values, [9,8,1,3])
        self.assertEqual(self.empty_col_end._coo_rows, [1,2,1,3])
        self.assertEqual(self.empty_col_end._coo_cols, [0,0,1,1])

        self.assertEqual(self.empty_col_end._pkd_ax, array([]))
        self.assertEqual(self.empty_col_end._unpkd_ax, array([]))
        self.assertEqual(self.empty_col_end._values, array([]))

        self.empty_cols.convert("csc")
        self.empty_cols.convert("coo")

        self.assertEqual(self.empty_cols._order, "coo")
        self.assertEqual(self.empty_cols._coo_values, [1,3])
        self.assertEqual(self.empty_cols._coo_rows, [1,3])
        self.assertEqual(self.empty_cols._coo_cols, [2,2])

        self.assertEqual(self.empty_cols._pkd_ax, array([]))
        self.assertEqual(self.empty_cols._unpkd_ax, array([]))
        self.assertEqual(self.empty_cols._values, array([]))

        # Test an empty matrix.
        self.empty.convert("csc")
        self.empty.convert("coo")

        self.assertEqual(self.empty._order, "coo")
        self.assertEqual(self.empty._coo_values, [])
        self.assertEqual(self.empty._coo_rows, [])
        self.assertEqual(self.empty._coo_cols, [])

        self.assertEqual(self.empty._pkd_ax, array([]))
        self.assertEqual(self.empty._unpkd_ax, array([]))
        self.assertEqual(self.empty._values, array([]))

    def test_convert_csc_csr(self):
        """convert csc to csr"""
        self.obj.convert("csc")
        self.obj.convert("csr")
        self.assertEqual(self.obj._order, "csr")
        self.assertEqual(self.obj._coo_values, [])
        self.assertEqual(self.obj._coo_rows, [])
        self.assertEqual(self.obj._coo_cols, [])

        self.assertEqual(self.obj._pkd_ax, [0, 2, 3, 4])
        self.assertEqual(self.obj._unpkd_ax, [0, 1, 2, 3])
        self.assertEqual(self.obj._values, [1,2,3,4])

        # Test a matrix with empty rows and columns.
        self.empty_cols.convert("csc")
        self.empty_cols.convert("csr")

        self.assertEqual(self.empty_cols._order, "csr")
        self.assertEqual(self.empty_cols._coo_values, [])
        self.assertEqual(self.empty_cols._coo_rows, [])
        self.assertEqual(self.empty_cols._coo_cols, [])

        self.assertEqual(self.empty_cols._pkd_ax, [0,0,1,1,2])
        self.assertEqual(self.empty_cols._unpkd_ax, [2,2])
        self.assertEqual(self.empty_cols._values, [1,3])

        # Test an empty matrix.
        self.empty.convert("csc")
        self.empty.convert("csr")

        self.assertEqual(self.empty._order, "csr")
        self.assertEqual(self.empty._coo_values, [])
        self.assertEqual(self.empty._coo_rows, [])
        self.assertEqual(self.empty._coo_cols, [])

        self.assertEqual(self.empty._pkd_ax, array([0,0,0,0]))
        self.assertEqual(self.empty._unpkd_ax, array([]))
        self.assertEqual(self.empty._values, array([]))

    def test_convert_csr_csc(self):
        """convert csr to csc"""
        self.obj.convert("csr")
        self.obj.convert("csc")
        self.assertEqual(self.obj._order, "csc")
        self.assertEqual(self.obj._coo_values, [])
        self.assertEqual(self.obj._coo_rows, [])
        self.assertEqual(self.obj._coo_cols, [])

        self.assertEqual(self.obj._pkd_ax, [0, 1, 2, 3, 4])
        self.assertEqual(self.obj._unpkd_ax, [0, 0, 1, 2])
        self.assertEqual(self.obj._values, [1,2,3,4])

        # Test a matrix with empty rows and columns.
        self.empty_cols.convert("csr")
        self.empty_cols.convert("csc")

        self.assertEqual(self.empty_cols._order, "csc")
        self.assertEqual(self.empty_cols._coo_values, [])
        self.assertEqual(self.empty_cols._coo_rows, [])
        self.assertEqual(self.empty_cols._coo_cols, [])

        self.assertEqual(self.empty_cols._pkd_ax, [0,0,0,2])
        self.assertEqual(self.empty_cols._unpkd_ax, [1,3])
        self.assertEqual(self.empty_cols._values, [1,3])

        # Test an empty matrix.
        self.empty.convert("csr")
        self.empty.convert("csc")

        self.assertEqual(self.empty._order, "csc")
        self.assertEqual(self.empty._coo_values, [])
        self.assertEqual(self.empty._coo_rows, [])
        self.assertEqual(self.empty._coo_cols, [])

        self.assertEqual(self.empty._pkd_ax, array([0,0,0,0,0]))
        self.assertEqual(self.empty._unpkd_ax, array([]))
        self.assertEqual(self.empty._values, array([]))

    def test_toCSR(self):
        """implicitly tested in test_convert_* functions"""
        pass

    def test_toCSC(self):
        """implicitly tested in test_convert_* functions"""
        pass

    def test_toCOO(self):
        """implicitly tested in test_convert_* functions"""
        pass

    def test_buildCSfromCOO(self):
        """implicitly tested in test_convert_* functions"""
        pass

    def test_buildCOOfromCS(self):
        """implicitly tested in test_convert_* functions"""
        pass

    def test_buildCSfromCS(self):
        """implicitly tested in test_convert_* functions"""
        pass

    def test_expand_compressed(self):
        """expand a compressed axis"""
        exp = [0,0,1,1,1,2,2,2,3,3]
        obs = self.obj._expand_compressed([0,2,5,8,10])
        self.assertEqual(obs,exp)

        exp = [1,1,1,1,1,2,2,2,3,3]
        obs = self.obj._expand_compressed([0,0,5,8,10])
        self.assertEqual(obs,exp)

        exp = [0,0,2,2,2,2,2,2,3,3]
        obs = self.obj._expand_compressed([0,2,2,8,10])
        self.assertEqual(obs,exp)

        exp = [0,0,1,1,1,2,2,2,2,2]
        obs = self.obj._expand_compressed([0,2,5,10,10])
        self.assertEqual(obs,exp)

        exp = [0,0,1,1,1,1,1,1,1,1]
        obs = self.obj._expand_compressed([0,2,10,10,10])
        self.assertEqual(obs,exp)

        exp = []
        obs = self.obj._expand_compressed([0,0,0,0,0])
        self.assertEqual(obs,exp)


class SupportTests(TestCase):
    def test_list_list_to_csmat(self):
        """convert [[row,col,value], ...] to csmat"""
        input = [[0,0,1],[1,1,5.0],[0,2,6]]
        exp = CSMat(2,3)
        exp.update({(0,0):1.0,(1,1):5.0,(0,2):6})
        obs = list_list_to_csmat(input)
        self.assertEqual(obs, exp)

    def test_nparray_to_csmat(self):
        """Convert nparray to csmat"""
        input = array([[1,2,3,4],[-1,6,7,8],[9,10,11,12]])
        exp = CSMat(3,4)
        exp.update({(0,0):1,(0,1):2,(0,2):3,(0,3):4,
                          (1,0):-1,(1,1):6,(1,2):7,(1,3):8,
                          (2,0):9,(2,1):10,(2,2):11,(2,3):12})
        obs = nparray_to_csmat(input)
        self.assertEqual(obs, exp)

    def test_list_dict_to_csmat(self):
        """Take a list of dicts and condense down to a single dict"""
        input = [{(0,0):10,(0,1):2}, {(1,2):15}, {(0,3):7}]
        exp = CSMat(3,4)
        exp.update({(0,0):10,(0,1):2,(1,2):15,(2,3):7})
        obs = list_dict_to_csmat(input)
        self.assertEqual(obs,exp)

    def test_dict_to_csmat(self):
        """Take a dict and convert to CSMat"""
        input = {(0,1):5,(1,0):2,(2,1):6}
        exp = CSMat(3,2)
        exp[(0,1)] = 5
        exp[(1,0)] = 2
        exp[(2,1)] = 6
        obs = dict_to_csmat(input)
        self.assertEqual(obs,exp)

    def test_to_csmat(self):
        """Convert to expected CSMat types"""
        vals = {(0,0):5,(0,1):6,(1,0):7,(1,1):8}
        obs = to_csmat(vals)
        exp = CSMat(2,2)
        exp[(0,0)] = 5
        exp[(0,1)] = 6
        exp[(1,0)] = 7
        exp[(1,1)] = 8
        self.assertEqual(obs,exp)

        input = {(0,1):5,(10,8):-1.23}

        exp = CSMat(11,9)
        exp[(0,1)] = 5
        exp[(10,8)] = -1.23
        obs = to_csmat(input)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # test transpose
        exp = CSMat(9,11)
        exp[(1,0)] = 5
        exp[(8,10)] = -1.23
        obs = to_csmat(input, transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of dicts, transpose
        exp = CSMat(3,2)
        exp[(0,0)] = 5.0
        exp[(1,0)] = 6.0
        exp[(2,0)] = 7.0
        exp[(0,1)] = 8.0
        exp[(1,1)] = 9.0
        exp[(2,1)] = 10.0
        obs = to_csmat([{(0,0):5,(0,1):6,(0,2):7},
                                           {(1,0):8,(1,1):9,(1,2):10}],
                                           transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of csmats
        exp = CSMat(2,3)
        exp[(0,0)] = 5
        exp[(0,1)] = 6
        exp[(0,2)] = 7
        exp[(1,0)] = 8
        exp[(1,1)] = 9
        exp[(1,2)] = 10
        row1 = CSMat(1,3)
        row1[(0,0)] = 5
        row1[(0,1)] = 6
        row1[(0,2)] = 7
        row2 = CSMat(1,3)
        row2[(0,0)] = 8
        row2[(0,1)] = 9
        row2[(0,2)] = 10
        obs = to_csmat([row1, row2])
        self.assertEqual(sorted(obs.items()), sorted(exp.items())) 

        # test empty set
        exp = CSMat(0,0)
        obs = to_csmat([])
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

    def test_list_nparray_to_csmat(self):
        """lists of nparrays to csmat"""
        ins = [array([0,2,1,0]), array([1,0,0,1])]
        exp = CSMat(2,4)
        exp[0,1] = 2
        exp[0,2] = 1
        exp[1,0] = 1
        exp[1,3] = 1
        obs = list_nparray_to_csmat(ins)
        self.assertEqual(obs,exp)

    def test_list_csmat_to_csmat(self):
        """list of csmats to csmat"""
        ins = [CSMat(1,4), CSMat(1,4)]
        ins[0][0,0] = 5
        ins[0][0,1] = 10
        ins[1][0,2] = 1
        ins[1][0,3] = 2
        exp = CSMat(2,4)
        exp[0,0] = 5
        exp[0,1] = 10
        exp[1,2] = 1
        exp[1,3] = 2
        obs = list_csmat_to_csmat(ins)
        self.assertEqual(obs,exp)


if __name__ == '__main__':
    main()
