#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from sparsemat import PySparseMatFloat, PySparseMatInt

class SparseMatFloatTests(TestCase):
    def setUp(self):
        self.obj = PySparseMatFloat()

    def test_get(self):
        """make sure we can get shibby"""
        self.obj.insert(10,20,30)
        self.assertEqual(self.obj.get(10,20), 30.0) # expect cast
        self.assertTrue(type(self.obj.get(10,20)), type(30.0))
        self.assertEqual(self.obj.get(10,19), 0.0)
        self.assertRaises(OverflowError, self.obj.get, -1, 10)
        self.assertRaises(OverflowError, self.obj.get, 1, -10)

    def test_insert(self):
        """make sure we can insert"""
        self.obj.insert(10,20,30)
        self.assertEqual(self.obj.get(10,20), 30.0) # expect cast
        self.assertTrue(type(self.obj.get(10,20)), type(30.0))
        self.obj.insert(10,20,-10.0)
        self.assertEqual(self.obj.get(10,20), -10.0)
        self.assertRaises(OverflowError, self.obj.insert, -1, 10, 2.0)
        self.assertRaises(OverflowError, self.obj.insert, 1, -10, 3.0)

    def test_getRow(self):
        self.fail()
    def test_getCol(self):
        self.fail()

    def test_eq(self):
        """equality..."""
        obj2 = PySparseMatFloat()
        self.assertEqual(self.obj, obj2)
        self.obj.insert(10,20,30.0)
        self.assertNotEqual(self.obj, obj2)

    def test_contains(self):
        self.fail()

if __name__ == '__main__':
    main()
