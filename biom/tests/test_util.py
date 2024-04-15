#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from os import remove
from os.path import abspath, dirname, exists
from tempfile import NamedTemporaryFile
from unittest import TestCase, main
import pathlib
import os

import numpy as np
import numpy.testing as npt

from biom.table import Table
from biom.parse import parse_biom_table, load_table
from biom.util import (natsort, flatten, unzip,
                       get_biom_project_dir, parse_biom_config_files,
                       compute_counts_per_sample_stats, safe_md5, biom_open,
                       get_data_path, generate_subsamples, is_hdf5_file)
import h5py

np.random.seed(1234)

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2017, The BIOM Format Development Team"
__credits__ = ["Rob Knight", "Peter Maxwell", "Sandra Smit",
               "Zongzhi Liu", "Micah Hamady", "Daniel McDonald",
               "Jai Ram Rideout", "Jorge Cañardo Alastuey"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class UtilTests(TestCase):

    def setUp(self):
        self.biom_otu_table1_w_tax = parse_biom_table(biom_otu_table1_w_tax)
        self.to_remove = []

    def tearDown(self):
        if self.to_remove:
            for f in self.to_remove:
                os.remove(f)

    def test_generate_subsamples(self):
        table = Table(np.array([[3, 1, 1], [0, 3, 3]]), ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual_o1 = set()
        actual_o2 = set()
        for i, obs in zip(range(100), generate_subsamples(table, 3)):
            actual_o1.add(tuple(obs.data('O1', 'observation')))
            actual_o2.add(tuple(obs.data('O2', 'observation')))

        self.assertEqual(actual_o1, {(3, 0, 0), (3, 1, 0), (3, 0, 1),
                                     (3, 1, 1)})
        self.assertEqual(actual_o2, {(0, 3, 3), (0, 2, 3), (0, 3, 2),
                                     (0, 2, 2)})

    def test_natsort(self):
        """natsort should perform numeric comparisons on strings

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        # string with alpha and numerics sort correctly
        s = 'sample1 sample2 sample11 sample12'.split()
        self.assertEqual(natsort(s),
                         'sample1 sample2 sample11 sample12'.split())
        s.reverse()
        self.assertEqual(natsort(s),
                         'sample1 sample2 sample11 sample12'.split())
        self.assertEqual(natsort(list('cba321')), list('123abc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdba')), list('abcd'))

        # string of ints sort correctly
        self.assertEqual(natsort(['11', '2', '1', '0']),
                         ['0', '1', '2', '11'])

        # strings of floats sort correctly
        self.assertEqual(natsort(['1.11', '1.12', '1.00', '0.009']),
                         ['0.009', '1.00', '1.11', '1.12'])

        # string of ints sort correctly
        self.assertEqual(
            natsort([('11', 'A'), ('2', 'B'), ('1', 'C'), ('0', 'D')]),
            [('0', 'D'), ('1', 'C'), ('2', 'B'), ('11', 'A')])

    def test_unzip(self):
        """unzip(items) should be the inverse of zip(*items)

        This method is ported from PyCogent (http://www.pycogent.org). PyCogent
        is a GPL project, but we obtained permission from the authors of this
        method to port it to the BIOM Format project (and keep it under BIOM's
        BSD license).
        """
        chars = [list('abcde'), list('ghijk')]
        numbers = [[1, 2, 3, 4, 5], [0, 0, 0, 0, 0]]
        strings = [["abcde", "fghij", "klmno"], ['xxxxx'] * 3]

        lists = [chars, numbers, strings]
        zipped = [zip(*i) for i in lists]
        unzipped = [unzip(i) for i in zipped]

        for u, l in zip(unzipped, lists):
            self.assertEqual(u, l)

    def test_flatten_no_change(self):
        """flatten should not change non-nested sequences (except to list)

        This method is ported from PyCogent (http://www.pycogent.org). PyCogent
        is a GPL project, but we obtained permission from the authors of this
        method to port it to the BIOM Format project (and keep it under BIOM's
        BSD license).
        """
        self.assertEqual(flatten('abcdef'), list('abcdef'))  # test identities
        self.assertEqual(flatten([]), [])  # test empty sequence
        self.assertEqual(flatten(''), [])  # test empty string

    def test_flatten(self):
        """flatten should remove one level of nesting from nested sequences

        This method is ported from PyCogent (http://www.pycogent.org). PyCogent
        is a GPL project, but we obtained permission from the authors of this
        method to port it to the BIOM Format project (and keep it under BIOM's
        BSD license).
        """
        self.assertEqual(flatten(['aa', 'bb', 'cc']), list('aabbcc'))
        self.assertEqual(flatten([1, [2, 3], [[4, [5]]]]), [1, 2, 3, [4, [5]]])

    def test_get_biom_project_dir(self):
        """Getting the biom project directory functions as expected.

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        # Do an explicit check on whether the file system containing
        # the current file is case insensitive. This is in response
        # to SF bug #2945548, where this test would fail on certain
        # unusual circumstances on case-insensitive file systems
        # because the case of abspath(__file__) was inconsistent.
        # (If you don't believe this, set case_insensitive_filesystem
        # to False, and rename your top-level biom-format directory as
        # Biom-format on OS X. That should cause this test to fail as
        # actual will be path/to/Biom-format and expected will be
        # path/to/biom-format.) Note that we don't need to change anything
        # in the get_biom_project_dir() function as if the
        # file system is case insenstive, the case of the returned
        # string is irrelevant.
        case_insensitive_filesystem = \
            exists(__file__.upper()) and exists(__file__.lower())

        actual = get_biom_project_dir()

        # I base the expected here off the imported location of
        # biom/util.py here, to handle cases where either the user
        # has biom-format in their PYTHONPATH, or when they've installed it
        # with setup.py.
        # If util.py moves this test will fail -- that
        # is what we want in this case, as the get_biom_project_dir()
        # function would need to be modified.
        import biom.util
        util_py_filepath = abspath(abspath(biom.util.__file__))
        expected = dirname(dirname(dirname(util_py_filepath)))

        if case_insensitive_filesystem:
            # Make both lowercase if the file system is case insensitive.
            actual = actual.lower()
            expected = expected.lower()
        self.assertEqual(actual, expected)

    def test_parse_biom_config_files(self):
        """parse_biom_config_files functions as expected.

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        fake_file1 = ['key1\tval1', 'key2 val2']
        fake_file2 = ['key2\tval3']
        actual = parse_biom_config_files([fake_file1, fake_file2])
        expected = {'key1': 'val1', 'key2': 'val3'}
        self.assertEqual(actual, expected)

        # Looking up a nonexistent value returns None.
        self.assertEqual(actual['fake_key'], None)

        # Empty dict on empty input.
        self.assertEqual(parse_biom_config_files([]), {})

    def test_compute_counts_per_sample_stats_empty(self):
        t = Table({}, [], [])
        res = compute_counts_per_sample_stats(t)
        self.assertEqual(res, (0, 0, 0, 0, {}))

    def test_compute_counts_per_sample_stats(self):
        """compute_counts_per_sample_stats functions as expected

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        actual = compute_counts_per_sample_stats(self.biom_otu_table1_w_tax)
        self.assertEqual(actual[0], 3)
        self.assertEqual(actual[1], 7)
        self.assertEqual(actual[2], 4)
        self.assertEqual(actual[3], 4.5)
        self.assertEqual(actual[4], {'Sample1': 7, 'Sample2': 3, 'Sample3': 4,
                                     'Sample4': 6, 'Sample5': 3, 'Sample6': 4})

    def test_compute_counts_per_sample_stats_obs_counts(self):
        """compute_counts_per_sample_stats functions as expected

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        actual = compute_counts_per_sample_stats(self.biom_otu_table1_w_tax,
                                                 binary_counts=True)
        self.assertEqual(actual[0], 1)
        self.assertEqual(actual[1], 4)
        self.assertEqual(actual[2], 2.5)
        self.assertEqual(actual[3], 2.5)
        self.assertEqual(actual[4], {'Sample1': 2, 'Sample2': 3, 'Sample3': 4,
                                     'Sample4': 2, 'Sample5': 1, 'Sample6': 3})

    def test_safe_md5(self):
        """Make sure we have the expected md5 with varied input types

        This method is ported from PyCogent (http://www.pycogent.org). PyCogent
        is a GPL project, but we obtained permission from the authors of this
        method to port it to the BIOM Format project (and keep it under BIOM's
        BSD license).
        """
        exp = 'd3b07384d113edec49eaa6238ad5ff00'

        tmp_f = NamedTemporaryFile(
            mode='w',
            prefix='test_safe_md5',
            suffix='txt',
            delete=False)
        tmp_f.write('foo\n')
        tmp_f.flush()

        obs = safe_md5(open(tmp_f.name))
        tmp_f.close()
        self.to_remove.append(tmp_f.name)
        self.assertEqual(obs, exp)

        obs = safe_md5(['foo\n'])
        self.assertEqual(obs, exp)

        # unsupported type raises TypeError
        self.assertRaises(TypeError, safe_md5, 42)

    def test_biom_open_hdf5_pathlib_write(self):
        t = Table(np.array([[0, 1, 2], [3, 4, 5]]), ['a', 'b'],
                  ['c', 'd', 'e'])
        with NamedTemporaryFile(delete=False) as tmpfile:
            with biom_open(pathlib.Path(tmpfile.name), 'w') as fp:
                t.to_hdf5(fp, 'tests')
        self.to_remove.append(tmpfile.name)

    def test_biom_open_hdf5_pathlib_read(self):
        cwd = os.getcwd()
        if '/' in __file__:
            os.chdir(__file__.rsplit('/', 1)[0])
        with biom_open(pathlib.Path('test_data/test.biom')) as f:
            self.assertTrue(isinstance(f, h5py.File))
        os.chdir(cwd)

    def test_biom_open_hdf5(self):
        with biom_open(get_data_path('test.biom')) as f:
            self.assertTrue(isinstance(f, h5py.File))

        with biom_open(get_data_path('test_writing.biom'), 'w') as f:
            self.assertTrue(isinstance(f, h5py.File))

        remove(get_data_path('test_writing.biom'))

    def test_biom_open_empty(self):
        with self.assertRaises(ValueError) as e:
            with biom_open(get_data_path('no-contents.biom'), 'r'):
                pass
        self.assertTrue("is empty and can't be parsed" in str(e.exception))

    def test_biom_open_json(self):
        with biom_open(get_data_path('test.json')) as f:
            self.assertTrue(hasattr(f, 'read'))

    def test_load_table_gzip_unicode(self):
        t = load_table(get_data_path('bad_table.txt.gz'))
        self.assertEqual('s__Cortinarius grosmornënsis',
                         t.metadata('otu1', 'observation')['taxonomy'])

    def test_load_table_unicode(self):
        t = load_table(get_data_path('bad_table.txt'))
        self.assertEqual('s__Cortinarius grosmornënsis',
                         t.metadata('otu1', 'observation')['taxonomy'])

    def test_is_hdf5_file(self):
        self.assertTrue(is_hdf5_file(get_data_path('test.biom')))
        self.assertFalse(is_hdf5_file(get_data_path('test.json')))

    def test_load_classic(self):
        tab = load_table(get_data_path('test.json'))
        with NamedTemporaryFile(mode='w', delete=False) as fp:
            fp.write(str(tab))
            fp.flush()

            obs = load_table(fp.name)
        self.to_remove.append(fp.name)

        npt.assert_equal(obs.ids(), tab.ids())
        npt.assert_equal(obs.ids(axis='observation'),
                         tab.ids(axis='observation'))
        npt.assert_equal(obs.matrix_data.toarray(), tab.matrix_data.toarray())


biom_otu_table1_w_tax = """{
     "id":null,
     "format": "Biological Observation Matrix 1.0.0-dev",
     "format_url": "http://biom-format.org",
     "type": "OTU table",
     "generated_by": "QIIME revision XYZ",
     "date": "2011-12-19T19:00:00",
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteoba\
cteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriac\
eae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_2", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobact\
eria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichosp\
ermum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchae\
ota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "\
g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicute\
s", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaer\
obium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}}
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
    }
"""


if __name__ == '__main__':
    main()
