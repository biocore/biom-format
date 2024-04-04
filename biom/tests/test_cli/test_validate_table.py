#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2017, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Daniel McDonald",
               "Jorge Cañardo Alastuey"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

import os
import json
from unittest import TestCase, main
from shutil import copy

import numpy as np

from biom.cli.table_validator import TableValidator
import h5py


class TableValidatorTests(TestCase):

    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = TableValidator()
        self.min_sparse_otu = json.loads(min_sparse_otu)
        self.rich_sparse_otu = json.loads(rich_sparse_otu)
        self.rich_dense_otu = json.loads(rich_dense_otu)
        self.min_dense_otu = json.loads(min_dense_otu)
        self.to_remove = []

        cur_path = os.path.split(os.path.abspath(__file__))[0]
        examples_path = os.path.join(cur_path.rsplit(os.path.sep, 3)[0],
                                     'examples')
        self.hdf5_file_valid = os.path.join(examples_path,
                                            'min_sparse_otu_table_hdf5.biom')
        self.hdf5_file_valid_md = os.path.join(examples_path,
                                               ('rich_sparse_otu_table_hdf5'
                                                '.biom'))

    def tearDown(self):
        for f in self.to_remove:
            os.remove(f)

    def test_valid_hdf5_metadata_v210(self):
        exp = {'valid_table': True, 'report_lines': []}
        obs = self.cmd(table=self.hdf5_file_valid,
                       format_version='2.1')
        self.assertEqual(obs, exp)
        obs = self.cmd(table=self.hdf5_file_valid_md,
                       format_version='2.1')
        self.assertEqual(obs, exp)

    def test_valid_hdf5_metadata_v200(self):
        pass  # omitting, not a direct way to test at this time using the repo

    def test_valid_hdf5(self):
        """Test a valid HDF5 table"""
        exp = {'valid_table': True,
               'report_lines': []}

        obs = self.cmd(table=self.hdf5_file_valid)
        self.assertEqual(obs, exp)

    def test_invalid_non_json(self):
        """Verify we error politely if a non-json ascii string is provided"""
        with self.assertRaisesRegex(ValueError,
                                    "^The provided table does not"):
            self.cmd(table=__file__)

    def test_invalid_hdf5(self):
        """Test an invalid HDF5 table"""
        exp = {'valid_table': False,
               'report_lines': ["Missing attribute: 'creation-date'"]}

        copy(self.hdf5_file_valid, 'invalid.hdf5')
        self.to_remove.append('invalid.hdf5')

        f = h5py.File('invalid.hdf5', 'a')
        del f.attrs['creation-date']

        f.close()
        obs = self.cmd(table='invalid.hdf5')
        self.assertEqual(obs, exp)

    def test_valid(self):
        """Correctly validates a table that is indeed... valid."""
        exp = {'valid_table': True, 'report_lines': []}

        f = open('valid_test1', 'w')
        f.write(json.dumps(self.min_sparse_otu))
        f.close()
        self.to_remove.append('valid_test1')

        obs = self.cmd(table='valid_test1')
        self.assertEqual(obs, exp)

        f = open('valid_test2', 'w')
        f.write(json.dumps(self.rich_sparse_otu))
        f.close()
        self.to_remove.append('valid_test2')

        obs = self.cmd(table='valid_test2')
        self.assertEqual(obs, exp)

        # Soldier, report!!
        f = open('valid_test3', 'w')
        f.write(json.dumps(self.rich_sparse_otu))
        f.close()
        self.to_remove.append('valid_test3')

        obs = self.cmd(table='valid_test3')
        self.assertTrue(obs['valid_table'])

    def test_invalid(self):
        """Correctly invalidates a table that is... invalid."""
        del self.min_sparse_otu['date']
        exp = {'valid_table': False, 'report_lines': ["Missing field: 'date'"]}

        f = open('invalid_test1', 'w')
        f.write(json.dumps(self.min_sparse_otu))
        f.close()
        self.to_remove.append('invalid_test1')

        obs = self.cmd(table='invalid_test1')
        self.assertEqual(obs, exp)

        self.rich_dense_otu['shape'][1] = 42
        exp = {'valid_table': False,
               'report_lines': ['Incorrect number of cols: [0, 0, 1, 0, 0, 0]',
                                "Number of columns in 'columns' is not equal "
                                "to 'shape'"]}

        f = open('invalid_test2', 'w')
        f.write(json.dumps(self.rich_dense_otu))
        f.close()
        self.to_remove.append('invalid_test2')

        obs = self.cmd(table='invalid_test2')
        self.assertEqual(obs, exp)

    def test_valid_format_url(self):
        """validates format url"""
        table = self.min_sparse_otu

        obs = self.cmd._valid_format_url(table)
        self.assertTrue(len(obs) == 0)

        table['format_url'] = 'foo'
        obs = self.cmd._valid_format_url(table)
        self.assertTrue(len(obs) > 0)

    def test_is_int(self):
        self.assertTrue(self.cmd._is_int(3))

        self.assertFalse(self.cmd._is_int(3.5))

        # checking with numpy dtypes
        self.assertFalse(self.cmd._is_int(np.float64(3)))
        self.assertFalse(self.cmd._is_int(np.float32(3)))

        self.assertTrue(self.cmd._is_int(np.int64(3)))
        self.assertTrue(self.cmd._is_int(np.int32(3)))
        self.assertTrue(self.cmd._is_int(np.int16(3)))

    def test_valid_format(self):
        """Should match format string"""
        table = self.min_sparse_otu

        self.cmd._format_version = '1.0.0'
        obs = self.cmd._valid_format(table)
        self.assertTrue(len(obs) == 0)

        table['format'] = 'foo'
        obs = self.cmd._valid_format(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_type(self):
        """Should be valid table type"""
        table = self.min_sparse_otu

        table['type'] = 'otu table'  # should not be case sensitive
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'Pathway table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'Function table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'Ortholog table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'Gene table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'Metabolite table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'OTU table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'Taxon table'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) == 0)

        table['type'] = 'foo'
        obs = self.cmd._valid_type(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_generated_by(self):
        """Should have some string for generated by"""
        table = self.min_sparse_otu
        obs = self.cmd._valid_generated_by(table)
        self.assertTrue(len(obs) == 0)

        table['generated_by'] = None
        obs = self.cmd._valid_generated_by(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_nullable_id(self):
        """Should just work."""
        pass

    def test_valid_metadata(self):
        """Can be nullable or an object"""
        table = self.min_sparse_otu

        table['rows'][2]['metadata'] = None
        obs = self.cmd._valid_metadata(table['rows'][2])
        self.assertTrue(len(obs) == 0)

        table['rows'][2]['metadata'] = {10: 20}
        obs = self.cmd._valid_metadata(table['rows'][2])
        self.assertTrue(len(obs) == 0)

        table['rows'][2]['metadata'] = ""
        obs = self.cmd._valid_metadata(table['rows'][2])
        self.assertTrue(len(obs) > 0)

        table['rows'][2]['metadata'] = "asdasda"
        obs = self.cmd._valid_metadata(table['rows'][2])
        self.assertTrue(len(obs) > 0)

        table['rows'][2]['metadata'] = [{'a': 'b'}, {'c': 'd'}]
        obs = self.cmd._valid_metadata(table['rows'][2])
        self.assertTrue(len(obs) > 0)

    def test_valid_matrix_type(self):
        """Make sure we have a valid matrix type"""
        obs = self.cmd._valid_matrix_type(self.min_dense_otu)
        self.assertTrue(len(obs) == 0)

        obs = self.cmd._valid_matrix_type(self.min_sparse_otu)
        self.assertTrue(len(obs) == 0)

        table = self.min_dense_otu

        table['matrix_type'] = 'spARSe'
        obs = self.cmd._valid_matrix_type(table)
        self.assertTrue(len(obs) > 0)

        table['matrix_type'] = 'sparse_asdasd'
        obs = self.cmd._valid_matrix_type(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_matrix_element_type(self):
        """Make sure we have a valid matrix type"""
        table = self.min_sparse_otu

        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_element_type'] = 'int'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_element_type'] = 'float'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_element_type'] = 'float'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_element_type'] = 'str'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_element_type'] = 'str'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_element_type'] = 'obj'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) > 0)

        table['matrix_element_type'] = 'asd'
        obs = self.cmd._valid_matrix_element_type(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_datetime(self):
        """Make sure we have a datetime stamp"""
        table = self.min_sparse_otu

        obs = self.cmd._valid_datetime(table)
        self.assertTrue(len(obs) == 0)

        table['date'] = "1999-11-11T10:11:12"
        obs = self.cmd._valid_datetime(table)
        self.assertTrue(len(obs) == 0)

    def test_valid_sparse_data(self):
        """Takes a sparse matrix field and validates"""
        table = self.min_sparse_otu

        obs = self.cmd._valid_sparse_data(table)
        self.assertTrue(len(obs) == 0)

        # incorrect type
        table['matrix_element_type'] = 'float'
        obs = self.cmd._valid_sparse_data(table)
        self.assertTrue(len(obs) > 0)

        # not balanced
        table['matrix_element_type'] = 'int'
        table['data'][5] = [0, 10]
        obs = self.cmd._valid_sparse_data(table)
        self.assertTrue(len(obs) > 0)

        # odd type for index
        table['data'][5] = [1.2, 5, 10]
        obs = self.cmd._valid_sparse_data(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_dense_data(self):
        """Takes a dense matrix field and validates"""
        table = self.min_dense_otu

        obs = self.cmd._valid_dense_data(table)
        self.assertTrue(len(obs) == 0)

        # incorrect type
        table['matrix_element_type'] = 'float'
        obs = self.cmd._valid_dense_data(table)
        self.assertTrue(len(obs) > 0)

        # not balanced
        table['matrix_element_type'] = 'int'
        table['data'][1] = [0, 10]
        obs = self.cmd._valid_dense_data(table)
        self.assertTrue(len(obs) > 0)

        # bad type in a field
        table['data'][1] = [5, 1, 0, 2.3, 3, 1]
        obs = self.cmd._valid_dense_data(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_shape(self):
        """validates shape information"""
        obs = self.cmd._valid_shape(self.min_sparse_otu)
        self.assertTrue(len(obs) == 0)

        obs = self.cmd._valid_shape(self.rich_sparse_otu)
        self.assertTrue(len(obs) == 0)

        bad_shape = self.min_sparse_otu.copy()
        bad_shape['shape'] = ['asd', 10]
        obs = self.cmd._valid_shape(bad_shape)
        self.assertTrue(len(obs) > 0)

    def test_valid_rows(self):
        """validates rows: field"""
        table = self.rich_dense_otu

        obs = self.cmd._valid_rows(table)
        self.assertTrue(len(obs) == 0)

        table['rows'][0]['id'] = ""
        obs = self.cmd._valid_rows(table)
        self.assertTrue(len(obs) > 0)

        table['rows'][0]['id'] = None
        obs = self.cmd._valid_rows(table)
        self.assertTrue(len(obs) > 0)

        del table['rows'][0]['id']
        obs = self.cmd._valid_rows(table)
        self.assertTrue(len(obs) > 0)

        table['rows'][0]['id'] = 'asd'
        table['rows'][0]['metadata'] = None
        obs = self.cmd._valid_rows(table)
        self.assertTrue(len(obs) == 0)

        # since this is an OTU table, metadata is a required key
        del table['rows'][0]['metadata']
        obs = self.cmd._valid_rows(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_columns(self):
        """validates table:columns: fields"""
        table = self.rich_dense_otu

        obs = self.cmd._valid_columns(table)
        self.assertTrue(len(obs) == 0)

        table['columns'][0]['id'] = ""
        obs = self.cmd._valid_columns(table)
        self.assertTrue(len(obs) > 0)

        table['columns'][0]['id'] = None
        obs = self.cmd._valid_columns(table)
        self.assertTrue(len(obs) > 0)

        del table['columns'][0]['id']
        obs = self.cmd._valid_columns(table)
        self.assertTrue(len(obs) > 0)

        table['columns'][0]['id'] = 'asd'
        table['columns'][0]['metadata'] = None
        obs = self.cmd._valid_columns(table)
        self.assertTrue(len(obs) == 0)

        # since this is an OTU table, metadata is a required key
        del table['columns'][0]['metadata']
        obs = self.cmd._valid_columns(table)
        self.assertTrue(len(obs) > 0)

    def test_valid_data(self):
        """validates data: fields"""
        # the burden of validating data is passed on to valid_sparse_data
        # and valid_dense_data
        table = self.rich_sparse_otu

        obs = self.cmd._valid_data(table)
        self.assertTrue(len(obs) == 0)

        table['matrix_type'] = 'foo'
        obs = self.cmd._valid_data(table)
        self.assertTrue(len(obs) > 0)


rich_sparse_otu = """{
"id":null,
"format": "1.0.0",
"format_url": "http://biom-format.org",
"type": "OTU table",
"generated_by": "QIIME revision XYZ",
"date": "2011-12-19T19:00:00",
"rows":[{"id":"GG_OTU_1",
         "metadata":{"taxonomy":["k__Bacteria",
                                 "p__Proteobacteria",
                                 "c__Gammaproteobacteria",
                                 "o__Enterobacteriales",
                                 "f__Enterobacteriaceae",
                                 "g__Escherichia",
                                 "s__"]}},
        {"id":"GG_OTU_2",
         "metadata":{"taxonomy":["k__Bacteria",
                                 "p__Cyanobacteria",
                                 "c__Nostocophycideae",
                                 "o__Nostocales",
                                 "f__Nostocaceae",
                                 "g__Dolichospermum",
                                 "s__"]}},
        {"id":"GG_OTU_3",
         "metadata":{"taxonomy":["k__Archaea",
                                 "p__Euryarchaeota",
                                 "c__Methanomicrobia",
                                 "o__Methanosarcinales",
                                 "f__Methanosarcinaceae",
                                 "g__Methanosarcina",
                                 "s__"]}},
        {"id":"GG_OTU_4",
         "metadata":{"taxonomy":["k__Bacteria",
                                 "p__Firmicutes",
                                 "c__Clostridia",
                                 "o__Halanaerobiales",
                                 "f__Halanaerobiaceae",
                                 "g__Halanaerobium",
                                 "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5",
         "metadata":{"taxonomy":["k__Bacteria",
                                 "p__Proteobacteria",
                                 "c__Gammaproteobacteria",
                                 "o__Enterobacteriales",
                                 "f__Enterobacteriaceae",
                                 "g__Escherichia",
                                 "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                             "BarcodeSequence":"CGCTTATCGAGA",
                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                             "BODY_SITE":"gut",
                             "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                             "BarcodeSequence":"CATACCAGTAGC",
                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                             "BODY_SITE":"gut",
                             "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                             "BarcodeSequence":"CTCTCTACCTGT",
                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                             "BODY_SITE":"gut",
                             "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                             "BarcodeSequence":"CTCTCGGCCTGT",
                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                             "BODY_SITE":"skin",
                             "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                             "BarcodeSequence":"CTCTCTACCAAT",
                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                             "BODY_SITE":"skin",
                             "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                             "BarcodeSequence":"CTAACTACCAAT",
                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                             "BODY_SITE":"skin",
                             "Description":"human skin"}}
        ],
     "matrix_type": "sparse",
     "matrix_element_type": "int",
     "shape": [5, 6],
     "data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,5,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
            ]
    }"""

min_sparse_otu = """{
        "id":null,
        "format": "1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "QIIME revision XYZ",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [5, 6],
        "data":[[0,2,1],
                [1,0,5],
                [1,1,1],
                [1,3,2],
                [1,4,3],
                [1,5,1],
                [2,2,1],
                [2,3,4],
                [2,5,2],
                [3,0,2],
                [3,1,1],
                [3,2,1],
                [3,5,1],
                [4,1,1],
                [4,2,1]
               ]
    }"""

rich_dense_otu = """{
     "id":null,
     "format": "1.0.0",
     "format_url": "http://biom-format.org",
     "type": "OTU table",
     "generated_by": "QIIME revision XYZ",
     "date": "2011-12-19T19:00:00",
     "rows":[{"id":"GG_OTU_1",
              "metadata":{"taxonomy":["k__Bacteria",
                                      "p__Proteobacteria",
                                      "c__Gammaproteobacteria",
                                      "o__Enterobacteriales",
                                      "f__Enterobacteriaceae",
                                      "g__Escherichia",
                                      "s__"]}},
             {"id":"GG_OTU_2",
              "metadata":{"taxonomy":["k__Bacteria",
                                      "p__Cyanobacteria",
                                      "c__Nostocophycideae",
                                      "o__Nostocales",
                                      "f__Nostocaceae",
                                      "g__Dolichospermum",
                                      "s__"]}},
             {"id":"GG_OTU_3",
              "metadata":{"taxonomy":["k__Archaea",
                                      "p__Euryarchaeota",
                                      "c__Methanomicrobia",
                                      "o__Methanosarcinales",
                                      "f__Methanosarcinaceae",
                                      "g__Methanosarcina",
                                      "s__"]}},
             {"id":"GG_OTU_4",
              "metadata":{"taxonomy":["k__Bacteria",
                                      "p__Firmicutes",
                                      "c__Clostridia",
                                      "o__Halanaerobiales",
                                      "f__Halanaerobiaceae",
                                      "g__Halanaerobium",
                                      "s__Halanaerobiumsaccharolyticum"]}},
             {"id":"GG_OTU_5",
              "metadata":{"taxonomy":["k__Bacteria",
                                      "p__Proteobacteria",
                                      "c__Gammaproteobacteria",
                                      "o__Enterobacteriales",
                                      "f__Enterobacteriaceae",
                                      "g__Escherichia",
                                      "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
                ],
     "matrix_type": "dense",
     "matrix_element_type": "int",
     "shape": [5,6],
     "data":  [[0,0,1,0,0,0],
               [5,1,0,2,3,1],
               [0,0,1,4,2,0],
               [2,1,1,0,0,1],
               [0,1,1,0,0,0]]
    }"""

min_dense_otu = """ {
        "id":null,
        "format": "1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "QIIME revision XYZ",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "dense",
        "matrix_element_type": "int",
        "shape": [5,6],
        "data":  [[0,0,1,0,0,0],
                  [5,1,0,2,3,1],
                  [0,0,1,4,2,0],
                  [2,1,1,0,0,1],
                  [0,1,1,0,0,0]]
    }"""


if __name__ == "__main__":
    main()
