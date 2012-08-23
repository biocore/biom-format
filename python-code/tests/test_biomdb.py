#!/usr/bin/env python
# encoding: utf-8

from cogent.util.unit_test import TestCase, main
from biomdb import BiomDB
from MySQLdb import ProgrammingError, OperationalError

class BiomDBTests(TestCase):
    def setUp(self):
        self.table_basename = 'biomdb_tests'
        self.biomdb = BiomDB(table_basename=self.table_basename)
    
    def tearDown(self):
        self.biomdb.dropTable()
        
    def test_tmp_table_count(self):
        """Count the number of tmp tables in the db"""
        exp = 0
        obs = self.biomdb._tmp_table_count()
        self.assertEqual(obs, exp)
        
        self.assertEqual(len(self.biomdb._tmp_tables), 0)
        
        newtable = self.biomdb.createTempTable()
        exp = 1 
        obs = self.biomdb._tmp_table_count()
        self.assertEqual(obs,exp)
        
        self.assertEqual(len(self.biomdb._tmp_tables), 1)
        
    def test_get_new_tableid(self):
        """Get a new table id"""
        exp = '_'.join([self.table_basename, '1'])
        obs = self.biomdb._get_new_tableid()
        self.assertEqual(obs,exp)
        
    def test_createTempTable_empty(self):
        """Create a new temporary table"""
        exp_name = '_'.join([self.table_basename, '1'])
        obs_name = self.biomdb.createTempTable()
        
        cursor = self.biomdb.con.cursor()
        cursor.execute("SELECT COUNT(*) FROM %s" % exp_name)
        exp = ((0,),)
        obs = cursor.fetchall()
        self.assertEqual(obs,exp)
        
    def test_createTempTable_fromexisting(self):
        """Creates a new temporary table from another table"""
        # create tmp table and put stuff in it
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, '1','a',10)
        self.biomdb.setItem(table, '10','c',20)
        
        new_table = self.biomdb.createTempTable(from_table=table)
        cursor = self.biomdb.con.cursor()
        cursor.execute("select * from %s" % new_table)
        obs = sorted(cursor.fetchall())
        exp = [('1','a',10), ('10','c',20)]
        self.assertEqual(obs,exp)
    
    def test_createTempTable_fromsamples(self):
        """Creates a new temporary table from another table samples"""
        # create tmp table and put stuff in it
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, '1','a',10)
        self.biomdb.setItem(table, '10','a',20)
        self.biomdb.setItem(table, '10','e',30)
        
        new_table = self.biomdb.createTempTable(from_table=table, 
                            from_samples=['a'])
        cursor = self.biomdb.con.cursor()
        cursor.execute("select * from %s" % new_table)
        obs = sorted(cursor.fetchall())
        exp = [('1','a',10),('10','a',20)]
        self.assertEqual(obs,exp)
        
        new_table = self.biomdb.createTempTable(table, ['x'])
        cursor = self.biomdb.con.cursor()
        cursor.execute("select * from %s" % new_table)
        obs = sorted(cursor.fetchall())
        exp = []
        self.assertEqual(obs,exp)
        
    def test_createTempTable_fromobs(self):
        """Creates a new temporary table from another table obs"""
        # create tmp table and put stuff in it
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, '1','a',10)
        self.biomdb.setItem(table, '10','a',20)
        self.biomdb.setItem(table, '10','e',30)
        
        new_table = self.biomdb.createTempTable(from_table=table, 
                            from_obs=['10'])
        cursor = self.biomdb.con.cursor()
        cursor.execute("select * from %s" % new_table)
        obs = sorted(cursor.fetchall())
        exp = [('10','a',20),('10','e',30)]
        self.assertEqual(obs,exp)
        
        new_table = self.biomdb.createTempTable(table, None, from_obs=['x'])
        cursor = self.biomdb.con.cursor()
        cursor.execute("select * from %s" % new_table)
        obs = sorted(cursor.fetchall())
        exp = []
        self.assertEqual(obs,exp)
        
    def test_createTempTable_bothaxes(self):
        """Bail if both axes are specified"""
        # create tmp table and put stuff in it
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, '1','a',10)
        self.biomdb.setItem(table, '10','a',20)
        self.biomdb.setItem(table, '10','e',30)
        
        self.assertRaises(ValueError, self.biomdb.createTempTable, table, ['a'],
                            ['10'])
                            
    def test_dropTable(self):
        """Drop temporary table or tables"""
        table = self.biomdb.createTempTable()
        self.biomdb.dropTable(table)
        self.assertRaises(ProgrammingError,self.biomdb.setItem,table,'a','b',10)
        self.assertRaises(OperationalError, self.biomdb.dropTable, table)
        
    def test_getItem(self):
        """gets an item"""
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, 'a','b',10)
        self.assertEqual(self.biomdb.getItem(table, 'a','b'), 10)
        self.assertEqual(self.biomdb.getItem(table, 'x','y'), 0)
        
    def test_setItem(self):
        """sets an item"""
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, 'a','b',10)
        self.assertEqual(self.biomdb.getItem(table, 'a','b'), 10)
        self.biomdb.setItem(table, 'a','b', 100)
        self.assertEqual(self.biomdb.getItem(table, 'a','b'), 100)
        cursor = self.biomdb.con.cursor()
        cursor.execute("select count(*) from %s" % table)
        exp = ((1L,),)
        obs = cursor.fetchall()
        self.assertEqual(obs,exp)
        
    def getSampleByObs(self):
        """Gets all samples associated to an observation"""
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, 'a','1',10)
        self.biomdb.setItem(table, 'a','2',20)
        self.biomdb.setItem(table, 'x','3',30)
        self.biomdb.setItem(table, 'x','4',40)
        
        exp = (('a','1',10), ('a','2',20))
        obs = sorted(self.biomdb.getSampleByObs(table, 'a'))
        self.assertEqual(obs,exp)
        
        exp = ()
        obs = self.biomdb.getSampleByObs(table,'asasd')
        self.assertEqual(obs,exp)
        
        
    def getObsBySample(self):
        """Gets all observations associated to a sample"""
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, 'a','1',10)
        self.biomdb.setItem(table, 'a','2',20)
        self.biomdb.setItem(table, 'x','2',30)
        self.biomdb.setItem(table, 'x','4',40)
        
        exp = (('a','2',20), ('x','2',30))
        obs = sorted(self.biomdb.getObsBySample(table, '2'))
        self.assertEqual(obs,exp)
        
        exp = ()
        obs = self.biomdb.getObsBySample(table,'asasd')
        self.assertEqual(obs,exp)
        
    def test_deleteItem(self):
        """deletes an item"""
        table = self.biomdb.createTempTable()
        self.biomdb.setItem(table, 'a','b',10)
        self.biomdb.deleteItem(table, 'a','b')
        self.assertEqual(self.biomdb.getItem(table,'a','b'), 0)
        
        # do it once more... dont complain if nothing to delete
        self.biomdb.deleteItem(table, 'a','b')
        self.assertEqual(self.biomdb.getItem(table,'a','b'), 0)
        
    def test_tableEquality_eq(self):
        """checks if tables are equal"""
        table1 = self.biomdb.createTempTable()
        self.biomdb.setItem(table1, '1','b',10)
        self.biomdb.setItem(table1, '2','b',20)
        self.biomdb.setItem(table1, '3','b',30)
        
        table2 = self.biomdb.createTempTable()
        self.biomdb.setItem(table2, '1','b',10)
        self.biomdb.setItem(table2, '2','b',20)
        self.biomdb.setItem(table2, '3','b',30)
        
        self.assertTrue(self.biomdb.tableEquality(table1, table2))
        
    def test_tableEquality_eq_value(self):
        """checks if tables are equal"""
        table1 = self.biomdb.createTempTable()
        self.biomdb.setItem(table1, '1','b',10)
        self.biomdb.setItem(table1, '2','b',20)
        self.biomdb.setItem(table1, '3','b',30)

        table2 = self.biomdb.createTempTable()
        self.biomdb.setItem(table2, '1','b',10)
        self.biomdb.setItem(table2, '2','b',20)
        self.biomdb.setItem(table2, '3','b',-30)

        self.assertFalse(self.biomdb.tableEquality(table1, table2))

    def test_tableEquality_eq_label(self):
        """checks if tables are equal"""
        table1 = self.biomdb.createTempTable()
        self.biomdb.setItem(table1, '1','b',10)
        self.biomdb.setItem(table1, '2','b',20)
        self.biomdb.setItem(table1, '3','b',30)

        table2 = self.biomdb.createTempTable()
        self.biomdb.setItem(table2, '1','b',10)
        self.biomdb.setItem(table2, '2','b',20)
        self.biomdb.setItem(table2, 'x','b',30)

        self.assertFalse(self.biomdb.tableEquality(table1, table2))
        
if __name__ == '__main__':
    main()