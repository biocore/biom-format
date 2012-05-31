#!/usr/bin/env python

from _sparsemat import PySparseMatFloat, PySparseMatInt
from numpy import zeros, ndarray, array
from biom.unit_test import TestCase, main
from biom.table import flatten
from biom.sparsemat import SparseMat, to_sparsemat, \
    list_nparray_to_sparsemat, list_list_to_sparsemat, \
    list_sparsemat_to_sparsemat, nparray_to_sparsemat, \
    dict_to_sparsemat, list_dict_to_sparsemat

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jai Rideout", "Greg Caporaso", 
               "Jose Clemente", "Justin Kuczynski"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "0.9.3-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
   
class SparseMatTests(TestCase):
    def setUp(self):
        self.obj = SparseMat(6,3)
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
        exp = SparseMat(1,3)
        exp[0,2] = 3
        self.assertEqual(exp, self.obj[1,:])

        exp = SparseMat(6,1)
        exp[1,0] = 3
        exp[5,0] = 6
        self.assertEqual(exp, self.obj[:,2])

        self.assertRaises(IndexError, self.obj.__getitem__, (10,slice(None)))
        self.assertRaises(AttributeError, self.obj.__getitem__, (3, slice(1,2,3)))

    def test_contains(self):
        """Make sure we can check things exist"""
        sm1 = SparseMat(3,4)
        for r in range(3):
            for c in range(4):
                assert (r,c) not in sm1
        sm1[1,2] = 0
        assert (1,2) not in sm1
        sm1[1,2] = 10
        assert (1,2) in sm1
        sm1.erase(1,2)
        assert (1,2) not in sm1
                
    def test_erase(self):
        """Make sure we can get rid of elements"""
        sm1 = SparseMat(4,6)
        self.assertEqual(sm1[2,3], 0.0)
        sm1.erase(2,3)
        self.assertEqual(sm1[2,3], 0.0)
        sm1[2,3] = 10
        self.assertEqual(sm1[2,3], 10.0)
        self.assertEqual(sm1._index_rows, [set([]), set([]), set([(2,3)]), set([])])
        self.assertEqual(sm1._index_cols, [set([]), set([]), set([]), set([(2,3)]), set([]), set([])])
        sm1.erase(2,3)
        self.assertEqual(sm1._index_rows, [set([]), set([]), set([]), set([])])
        self.assertEqual(sm1._index_cols, [set([]), set([]), set([]), set([]), set([]), set([])])
        self.assertEqual(sm1[2,3], 0.0)
        self.assertEqual(sm1._index_rows, [set([]), set([]), set([]), set([])])
        self.assertEqual(sm1._index_cols, [set([]), set([]), set([]), set([]), set([]), set([])])
        
    def test_eq(self):
        """Tests for equality"""
        sm1 = SparseMat(4,6)
        sm2 = SparseMat(4,6)
        sm3 = SparseMat(6,4)
        
        self.assertEqual(sm1, sm2)
        self.assertNotEqual(sm1,sm3)
        
        sm1[0,1] = 10
        sm2[0,1] = 5
        self.assertNotEqual(sm1, sm2)
        sm2[0,1] = 10
        self.assertEqual(sm1, sm2)
        
    def test_update_internal_indices(self):
        """Update internal indices"""
        sd = SparseMat(2,3)
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
        exp = SparseMat(1,3)
        exp.update({(0,2):3})
        obs = self.obj.getRow(1)
        self.assertEqual(obs, exp)

        exp = SparseMat(1,3)
        obs = self.obj.getRow(4)
        self.assertEqual(obs,exp)

        self.assertRaises(IndexError, self.obj.getRow, -1)

    def test_getCol(self):
        """Get a col"""
        exp = SparseMat(6,1)
        exp.update({(1,0):3,(5,0):6})
        obs = self.obj.getCol(2)
        self.assertEqual(obs,exp)

        exp = SparseMat(6,1)
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
        
        self.assertEqual(sorted(self.obj.items()), sorted({(1,2):1.0,(1,1):10.0,(5,2):6.0}.items()))
        self.assertEqual(self.obj._index_rows, [set(),set([(1,2),(1,1)]),\
                                                set(),set(),\
                                                set(), set([(5,2)])])
        self.assertEqual(self.obj._index_cols, [set(),set([(1,1)]),set([(1,2),(5,2)])])

    def test_T(self):
        """test transpose"""
        exp = SparseMat(3,6)
        exp.update({(2,1):3,(2,5):6})
        obs = self.obj.T
        self.assertEqual(obs, exp)

class SupportTests(TestCase):
    def test_list_list_to_sparsemat(self):
        """convert [[row,col,value], ...] to sparsemat"""
        input = [[0,0,1],[10,0,5.0],[2,3,6]]
        exp = SparseMat(11,4)
        exp.update({(0,0):1.0,(10,0):5.0,(2,3):6})
        obs = list_list_to_sparsemat(input)
        self.assertEqual(obs, exp)

    def test_nparray_to_sparsemat(self):
        """Convert nparray to sparsemat"""
        input = array([[1,2,3,4],[-1,6,7,8],[9,10,11,12]])
        exp = SparseMat(3,4)
        exp.update({(0,0):1,(0,1):2,(0,2):3,(0,3):4,
                          (1,0):-1,(1,1):6,(1,2):7,(1,3):8,
                          (2,0):9,(2,1):10,(2,2):11,(2,3):12})
        obs = nparray_to_sparsemat(input)
        self.assertEqual(obs, exp)

    def test_list_dict_to_sparsemat(self):
        """Take a list of dicts and condense down to a single dict"""
        input = [{(0,5):10,(10,10):2}, {(0,1):15}, {(0,3):7}]
        exp = SparseMat(3,11)
        exp.update({(0,5):10,(0,10):2,(1,1):15,(2,3):7})
        obs = list_dict_to_sparsemat(input)
        self.assertEqual(obs,exp)

    def test_dict_to_sparsemat(self):
        """Take a dict and convert to SparseMat"""
        input = {(0,10):5,(4,5):2,(3,3):6}
        exp = SparseMat(5,11)
        exp[(0,10)] = 5
        exp[(4,5)] = 2
        exp[(3,3)] = 6
        obs = dict_to_sparsemat(input)
        self.assertEqual(obs,exp)

    def test_to_sparsemat(self):
        """Convert to expected SparseMat types"""
        vals = {(0,0):5,(0,1):6,(1,0):7,(1,1):8}
        obs = to_sparsemat(vals)
        exp = SparseMat(2,2)
        exp[(0,0)] = 5
        exp[(0,1)] = 6
        exp[(1,0)] = 7
        exp[(1,1)] = 8
        self.assertEqual(obs,exp)

        input = {(0,1):5,(10,8):-1.23}

        exp = SparseMat(11,9)
        exp[(0,1)] = 5
        exp[(10,8)] = -1.23
        obs = to_sparsemat(input)
        self.assertEqual(obs.items(), exp.items())

        # test transpose
        exp = SparseMat(9,11)
        exp[(1,0)] = 5
        exp[(8,10)] = -1.23
        obs = to_sparsemat(input, transpose=True)
        self.assertEqual(obs.items(), exp.items())

        # passing a list of dicts, transpose
        exp = SparseMat(3,2)
        exp[(0,0)] = 5.0
        exp[(1,0)] = 6.0
        exp[(2,0)] = 7.0
        exp[(0,1)] = 8.0
        exp[(1,1)] = 9.0
        exp[(2,1)] = 10.0
        obs = to_sparsemat([{(0,0):5,(0,1):6,(0,2):7},
                                           {(1,0):8,(1,1):9,(1,2):10}],
                                           transpose=True)
        self.assertEqual(sorted(obs.items()), sorted(exp.items()))

        # passing a list of sparsemats
        exp = SparseMat(2,3)
        exp[(0,0)] = 5
        exp[(0,1)] = 6
        exp[(0,2)] = 7
        exp[(1,0)] = 8
        exp[(1,1)] = 9
        exp[(1,2)] = 10
        row1 = SparseMat(1,3)
        row1[(0,0)] = 5
        row1[(0,1)] = 6
        row1[(0,2)] = 7
        row2 = SparseMat(1,3)
        row2[(0,0)] = 8
        row2[(0,1)] = 9
        row2[(0,2)] = 10
        obs = to_sparsemat([row1, row2])
        self.assertEqual(obs.items(), exp.items()) 

        # test empty set
        exp = SparseMat(0,0)
        obs = to_sparsemat([])
        self.assertEqual(obs.items(), exp.items())

class PySparseMatIntTests(TestCase):
    def setUp(self):
        self.obj = PySparseMatInt(20,30)

    def test_get(self):
        """make sure we can get shibby"""
        self.obj.insert(10,20,30)
        self.assertEqual(self.obj.get(10,20), 30) # expect cast
        self.assertTrue(type(self.obj.get(10,20)), type(30))
        self.assertEqual(self.obj.get(10,19), 0)
        self.assertRaises(KeyError, self.obj.get, -1, 10)
        self.assertRaises(KeyError, self.obj.get, 1, -10)

    def test_insert(self):
        """make sure we can insert"""
        self.obj.insert(10,20,30)
        self.assertEqual(self.obj.get(10,20), 30) # expect cast
        self.assertTrue(type(self.obj.get(10,20)), type(30))
        self.obj.insert(10,20,-10)
        self.assertEqual(self.obj.get(10,20), -10)
        self.assertRaises(KeyError, self.obj.insert, -1, 10, 2)
        self.assertRaises(KeyError, self.obj.insert, 1, -10, 3)

    def test_contains(self):
        """Make sure we can check if things are present"""
        x = PySparseMatInt(2,3)
        self.assertEqual(x.contains(1,2), 0)
        x.insert(1,2,20)
        self.assertEqual(x.contains(1,2), 1)
        x.insert(1,2,10)
        self.assertEqual(x.contains(1,2), 1)
        x.erase(1,2)
        self.assertEqual(x.contains(1,2), 0)
        x.insert(1,2,0.0)
        self.assertEqual(x.contains(1,2), 0)
        
    def test_length(self):
        """make sure we can test length"""
        x = PySparseMatInt(3,4)
        self.assertEqual(x.length(), 0)
        x.insert(1,2,10)
        self.assertEqual(x.length(), 1)
        x.insert(2,3,4)
        x.insert(2,3,4)
        x.insert(2,3,4)
        self.assertEqual(x.length(), 2)
        x.erase(2,3)
        self.assertEqual(x.length(), 1)
        
    def test_erase(self):
        """make sure we can erase"""
        x = PySparseMatInt(2,3)
        x.insert(1,2,10)
        self.assertEqual(x.get(1,2), 10)
        x.erase(1,2)
        self.assertEqual(x.get(1,2), 0)
        self.assertEqual(x.contains(1,2), 0)
    
    def test_keys(self):
        """make sure we can get keys"""
        x = PySparseMatInt(3,4)
        self.assertEqual(x.keys(), [])
        x.insert(1,2,10)
        x.insert(2,3,4)
        x.insert(2,3,4)
        x.insert(2,3,4)
        self.assertEqual(sorted(x.keys()), [(1,2),(2,3)])
        x.erase(2,3)
        self.assertEqual(x.keys(), [(1,2)])
            
    def test_items(self):
        """make sure we can get items"""
        x = PySparseMatInt(3,4)
        self.assertEqual(x.items(), [])
        x.insert(1,2,10)
        x.insert(2,3,4)
        x.insert(2,3,4)
        x.insert(2,3,4)
        self.assertEqual(sorted(x.items()), [((1,2),10),((2,3),4)])
        x.erase(2,3)
        self.assertEqual(x.items(), [((1,2),10)])
        
class PySparseMatFloatTests(TestCase):
    def setUp(self):
        self.obj = PySparseMatFloat(11,21)

    def test_get(self):
        """make sure we can get shibby"""
        self.obj.insert(10,20,30)
        self.assertEqual(self.obj.get(10,20), 30.0) # expect cast
        self.assertTrue(type(self.obj.get(10,20)), type(30.0))
        self.assertEqual(self.obj.get(10,19), 0.0)
        self.assertRaises(KeyError, self.obj.get, -1, 10)
        self.assertRaises(KeyError, self.obj.get, 1, -10)

    def test_insert(self):
        """make sure we can insert"""
        self.obj.insert(10,20,30)
        self.assertEqual(self.obj.get(10,20), 30.0) # expect cast
        self.assertTrue(type(self.obj.get(10,20)), type(30.0))
        self.obj.insert(10,20,-10.0)
        self.assertEqual(self.obj.get(10,20), -10.0)
        self.assertRaises(KeyError, self.obj.insert, -1, 10, 2.0)
        self.assertRaises(KeyError, self.obj.insert, 1, -10, 3.0)

    def test_contains(self):
        """Make sure we can check if things are present"""
        x = PySparseMatFloat(2,3)
        self.assertEqual(x.contains(1,2), 0)
        x.insert(1,2,20)
        self.assertEqual(x.contains(1,2), 1)
        x.insert(1,2,10)
        self.assertEqual(x.contains(1,2), 1)
        x.erase(1,2)
        self.assertEqual(x.contains(1,2), 0)
        x.insert(1,2,0.0)
        self.assertEqual(x.contains(1,2), 0)
        
    def test_length(self):
        """make sure we can test length"""
        x = PySparseMatFloat(3,4)
        self.assertEqual(x.length(), 0)
        x.insert(1,2,10)
        self.assertEqual(x.length(), 1)
        x.insert(2,3,4)
        x.insert(2,3,4)
        x.insert(2,3,4)
        self.assertEqual(x.length(), 2)
        x.erase(2,3)
        self.assertEqual(x.length(), 1)
        
    def test_erase(self):
        """make sure we can erase"""
        x = PySparseMatFloat(2,3)
        x.insert(1,2,10)
        self.assertEqual(x.get(1,2), 10.0)
        x.erase(1,2)
        self.assertEqual(x.get(1,2), 0.0)
        self.assertEqual(x.contains(1,2), 0)
    
    def test_keys(self):
        """make sure we can get keys"""
        x = PySparseMatFloat(3,4)
        self.assertEqual(x.keys(), [])
        x.insert(1,2,10)
        x.insert(2,3,4)
        x.insert(2,3,4)
        x.insert(2,3,4)
        self.assertEqual(sorted(x.keys()), [(1,2),(2,3)])
        x.erase(2,3)
        self.assertEqual(x.keys(), [(1,2)])
            
    def test_items(self):
        """make sure we can get items"""
        x = PySparseMatFloat(3,4)
        self.assertEqual(x.items(), [])
        x.insert(1,2,10)
        x.insert(2,3,4)
        x.insert(2,3,4)
        x.insert(2,3,4)
        self.assertEqual(sorted(x.items()), [((1,2),10.0),((2,3),4.0)])
        x.erase(2,3)
        self.assertEqual(x.items(), [((1,2),10.0)])
        
if __name__ == '__main__':
    main()
