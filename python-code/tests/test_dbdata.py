#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from dbdata import dbData, biomdb
from numpy import array

class dbdataTests(TestCase):
    def setUp(self):
        self.basic = dbData(['a','b','c','d'],['1','2','3','4'])
    
    def test_init(self):
        """make sure shibby is constructed all proper like"""
        self.assertEqual(self.basic.shape, (4,4))
        self.assertEqual(self.basic.dtype, 'FLOAT')
        self.assertEqual(self.basic.Table.rsplit('_',1)[0], 'tmp_biomtable') 
        
        o_r_i = {0:'a',1:'b',2:'c',3:'d'}
        o_i = {'a':0,'b':1,'c':2,'d':3}
        s_r_i = {0:'1',1:'2',2:'3',3:'4'}
        s_i = {'1':0,'2':1,'3':2,'4':3}
        self.assertEqual(self.basic._obs_r_index, o_r_i)
        self.assertEqual(self.basic._obs_index, o_i)
        self.assertEqual(self.basic._sample_r_index, s_r_i)
        self.assertEqual(self.basic._sample_index, s_i)
        
    def test_setitem(self):
        """it puts the items into the db"""
        self.basic[0,0] = 100
        self.basic[0,2] = 200.0
        self.assertEqual(self.basic[0,0], 100.0)
        self.assertEqual(self.basic[0,2], 200.0)
        self.basic[0,0] = 0
        self.assertEqual(self.basic[0,0], 0)
        self.assertRaises(IndexError, self.basic.__setitem__, (-1,0), 2)
        self.assertRaises(IndexError, self.basic.__setitem__, (0,-1), 2)
        self.assertRaises(IndexError, self.basic.__setitem__, (4,3), 2)
        self.assertRaises(IndexError, self.basic.__setitem__, (3,4), 2)
        
    def test_getitem(self):
        """it gets the items"""
        self.basic[0,0] = 100
        self.basic[0,2] = 200.0
        self.assertEqual(self.basic[0,0], 100.0)
        self.assertEqual(self.basic[0,2], 200.0)
        self.basic[0,0] = 0
        self.assertEqual(self.basic[0,0], 0)
        self.assertRaises(IndexError, self.basic.__getitem__, (-1,0))
        self.assertRaises(IndexError, self.basic.__getitem__, (0,-1))
        self.assertRaises(IndexError, self.basic.__getitem__, (4,3))
        self.assertRaises(IndexError, self.basic.__getitem__, (3,4))
        
        # slicing wraps getRow/getCol, tests in there
        
    def test_getRow(self):
        """it gets all the things in a row"""
        self.basic[0,0] = 100
        self.basic[0,2] = 200.0
        exp = array([100.0,0,200.0,0])
        obs = self.basic.getRow(0)
        
        exp = array([0,0,0,0])
        obs = self.basic.getRow(1)
        
        self.assertRaises(IndexError, self.basic.getRow, -1)
        self.assertRaises(IndexError, self.basic.getRow, 5)
        
    def test_getCol(self):
        """it gets all the things in a col"""
        self.basic[0,0] = 100
        self.basic[2,0] = 200.0
        exp = array([100.0,0,200.0,0])
        obs = self.basic.getCol(0)
        
        exp = array([0,0,0,0])
        obs = self.basic.getCol(1)
        
        self.assertRaises(IndexError, self.basic.getCol, -1)
        self.assertRaises(IndexError, self.basic.getCol, 5)
    
    def test_copy(self):
        """copy all the things"""
        self.basic[0,0] = 100
        self.basic[1,2] = 2.32
        new_table = self.basic.copy()
        
        self.assertNotEqual(self.basic.Table, new_table.Table)
        self.assertEqual(self.basic.dtype, new_table.dtype)
        self.assertEqual(self.basic, new_table)
        self.assertNotSameObj(self.basic, new_table)
        self.assertEqual(new_table[0,0], 100)
        self.assertEqual(new_table[1,2], 2.32)
        
    def test_transpose(self):
        """transpose all the things"""
        table = dbData(['1','2','3'],['a','b'])
        table[0,0] = 5
        table[0,1] = 10
        table[2,1] = 15
        
        print "before transpose"
        transposed = table.transpose()
        print "after transpose"
        self.assertEqual(transposed.shape, (2,3))
        self.assertEqual(transposed.Observations, ['a','b'])
        self.assertEqual(transposed.Samples, ['1','2','3'])
        self.assertEqual(transposed[0,0], 5)
        self.assertEqual(transposed[1,0], 10)
        self.assertEqual(transposed[1,2], 15)

if __name__ == '__main__':
    main()