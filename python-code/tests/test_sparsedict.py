#!/usr/bin/env python

from numpy import zeros, ndarray, array
from biom.unit_test import TestCase, main
from biom.table import flatten
from biom.sparsedict import SparseDict, to_sparsedict, \
    list_nparray_to_sparsedict, list_list_to_sparsedict, \
    list_dict_to_sparsedict, nparray_to_sparsedict, \
    dict_to_sparsedict
 
__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso", 
               "Jose Clemente", "Justin Kuczynski"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
    
class SparseDictTests(TestCase):
    def setUp(self):
        self.obj = SparseDict(6,3)
        self.obj.update({(1,2):3,(5,2):6})

    def test_copy(self):
        """copy thy self"""
        obs = self.obj.copy()
        self.assertEqual(obs, self.obj)
        self.assertNotEqual(id(obs), id(self.obj))
        obs[1,2] = 10
        self.assertNotEqual(obs, self.obj)
        obs[1,2] = 3
        obs[0,0] = 5
        self.assertNotEqual(obs, self.obj)
        self.assertNotEqual(id(self.obj._index_rows), id(obs._index_rows))
        self.assertNotEqual(id(self.obj._index_cols), id(obs._index_cols))

    def test_setitem(self):
        self.obj[(2,2)] = 10
        exp = sorted([((1,2),3),((5,2),6),((2,2),10)])
        self.assertEqual(sorted(self.obj.items()), exp)
        self.assertRaises(KeyError, self.obj.__setitem__, (100,50), 10)

        self.assertEqual(self.obj._index_rows, [set(),set([(1,2)]),
                                                set([(2,2)]), set(), set(),
                                                set([(5,2)])])
        self.assertEqual(self.obj._index_cols, [set(),set(),
                                                set([(1,2),(2,2),(5,2)])])


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
        exp = SparseDict(1,3)
        exp[0,2] = 3
        self.assertEqual(exp, self.obj[1,:])

        exp = SparseDict(6,1)
        exp[1,0] = 3
        exp[5,0] = 6
        self.assertEqual(exp, self.obj[:,2])

        self.assertRaises(IndexError, self.obj.__getitem__, (10,slice(None)))
        self.assertRaises(AttributeError, self.obj.__getitem__, (3, slice(1,2,3)))

    def test_update_internal_indices(self):
        """Update internal indices"""
        sd = SparseDict(2,3)
        self.assertEqual(sd._index_rows, [set(),set()])
        self.assertEqual(sd._index_cols, [set(),set(),set()])

        sd[(1,2)] = 5
        self.assertEqual(sd._index_rows, [set(),set([(1,2)])])
        self.assertEqual(sd._index_cols, [set(),set(),set([(1,2)])])

        sd[(1,2)] = 0
        self.assertEqual(sd._index_rows, [set(),set()])
        self.assertEqual(sd._index_cols, [set(),set(),set()])

        sd[(1,1)] = 0
        self.assertEqual(sd._index_rows, [set(),set()])
        self.assertEqual(sd._index_cols, [set(),set(),set()])

    def test_getRow(self):
        """Get a row"""
        exp = SparseDict(1,3)
        exp.update({(0,2):3})
        obs = self.obj.getRow(1)
        self.assertEqual(obs, exp)

        exp = SparseDict(1,3)
        obs = self.obj.getRow(4)
        self.assertEqual(obs,exp)

        self.assertRaises(IndexError, self.obj.getRow, -1)

    def test_getCol(self):
        """Get a col"""
        exp = SparseDict(6,1)
        exp.update({(1,0):3,(5,0):6})
        obs = self.obj.getCol(2)
        self.assertEqual(obs,exp)

        exp = SparseDict(6,1)
        obs = self.obj.getCol(1)
        self.assertEqual(obs,exp)

        self.assertRaises(IndexError, self.obj.getCol, -1)

    def test_update(self):
        """updates should work and update indexes"""
        items = self.obj.items()
        indexes = (self.obj._index_rows, self.obj._index_cols)
        self.obj.update({(1,2):3,(5,2):6})
        self.assertEqual(items, self.obj.items())
        self.assertEqual(indexes, (self.obj._index_rows, self.obj._index_cols))

        self.obj.update({(1,2):0,(5,2):6})
        self.assertEqual(self.obj.items(), {(5,2):6}.items())
        self.assertEqual(self.obj._index_rows, [set(),set(),set(),set(),\
                                                set(), set([(5,2)])])
        self.assertEqual(self.obj._index_cols, [set(),set(),set([(5,2)])])

        self.obj.update({(1,2):1,(2,2):0,(1,1):10})
        self.assertEqual(self.obj.items(), {(1,2):1,(1,1):10,(5,2):6}.items())
        self.assertEqual(self.obj._index_rows, [set(),set([(1,2),(1,1)]),\
                                                set(),set(),\
                                                set(), set([(5,2)])])
        self.assertEqual(self.obj._index_cols, [set(),set([(1,1)]),set([(1,2),(5,2)])])

    def test_T(self):
        """test transpose"""
        exp = SparseDict(3,6)
        exp.update({(2,1):3,(2,5):6})
        obs = self.obj.T
        self.assertEqual(obs, exp)

    def test_get_size(self):
        """test getting the number of nonzero elements"""
        self.assertEqual(self.obj.size, 2)

        # Test with setting an element explicitly to zero.
        sd = SparseDict(2,4)
        sd.update({(0,1):3,(1,2):7,(0,0):0})
        self.assertEqual(sd.size, 2)

        # Test with an empty matrix.
        sd = SparseDict(2,4)
        self.assertEqual(sd.size, 0)
        sd = SparseDict(0,0)
        self.assertEqual(sd.size, 0)


class SupportTests(TestCase):
    def test_list_list_to_sparsedict(self):
        """convert [[row,col,value], ...] to dict"""
        input = [[0,0,1],[10,0,5.0],[2,3,6]]
        exp = {(0,0):1.0,(10,0):5.0,(2,3):6}
        obs = list_list_to_sparsedict(input)
        self.assertEqual(obs, exp)

    def test_nparray_to_sparsedict(self):
        """Convert nparray to dict"""
        input = array([[1,2,3,4],[-1,6,7,8],[9,10,11,12]])
        exp = SparseDict(3,4)
        exp.update({(0,0):1,(0,1):2,(0,2):3,(0,3):4,
                          (1,0):-1,(1,1):6,(1,2):7,(1,3):8,
                          (2,0):9,(2,1):10,(2,2):11,(2,3):12})
        obs = nparray_to_sparsedict(input)
        self.assertEqual(obs, exp)

    def test_list_dict_to_sparsedict(self):
        """Take a list of dicts and condense down to a single dict"""
        input = [{(0,5):10,(10,10):2}, {(0,1):15}, {(0,3):7}]
        exp = SparseDict(3,11)
        exp.update({(0,5):10,(0,10):2,(1,1):15,(2,3):7})
        obs = list_dict_to_sparsedict(input)
        self.assertEqual(obs,exp)

    def test_dict_to_sparsedict(self):
        """Take a dict and convert to SparseDict"""
        input = {(0,10):5,(4,5):2,(3,3):6}
        exp = SparseDict(5,11)
        exp[(0,10)] = 5
        exp[(4,5)] = 2
        exp[(3,3)] = 6
        obs = dict_to_sparsedict(input)
        self.assertEqual(obs,exp)

    def test_to_sparsedict(self):
        """Convert to expected SparseDict types"""
        vals = {(0,0):5,(0,1):6,(1,0):7,(1,1):8}
        obs = to_sparsedict(vals)
        exp = SparseDict(2,2)
        exp[(0,0)] = 5
        exp[(0,1)] = 6
        exp[(1,0)] = 7
        exp[(1,1)] = 8
        self.assertEqual(obs,exp)

        input = {(0,1):5,(10,8):-1.23}

        exp = SparseDict(11,9)
        exp[(0,1)] = 5
        exp[(10,8)] = -1.23
        obs = to_sparsedict(input)
        self.assertEqual(obs.items(), exp.items())

        # test transpose
        exp = SparseDict(9,11)
        exp[(1,0)] = 5
        exp[(8,10)] = -1.23
        obs = to_sparsedict(input, transpose=True)
        self.assertEqual(obs.items(), exp.items())

        # passing a list of dicts, transpose
        exp = SparseDict(3,2)
        exp[(0,0)] = 5.0
        exp[(1,0)] = 6.0
        exp[(2,0)] = 7.0
        exp[(0,1)] = 8.0
        exp[(1,1)] = 9.0
        exp[(2,1)] = 10.0
        obs = to_sparsedict([{(0,0):5,(0,1):6,(0,2):7},
                                           {(1,0):8,(1,1):9,(1,2):10}],
                                           transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of sparsedicts
        exp = SparseDict(2,3)
        exp[(0,0)] = 5
        exp[(0,1)] = 6
        exp[(0,2)] = 7
        exp[(1,0)] = 8
        exp[(1,1)] = 9
        exp[(1,2)] = 10
        row1 = SparseDict(1,3)
        row1[(0,0)] = 5
        row1[(0,1)] = 6
        row1[(0,2)] = 7
        row2 = SparseDict(1,3)
        row2[(0,0)] = 8
        row2[(0,1)] = 9
        row2[(0,2)] = 10
        obs = to_sparsedict([row1, row2])
        self.assertEqual(obs.items(), exp.items()) 

        # test empty set
        exp = SparseDict(0,0)
        obs = to_sparsedict([])
        self.assertEqual(obs.items(), exp.items())

if __name__ == '__main__':
    main()
