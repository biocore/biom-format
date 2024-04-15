#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import os
from io import StringIO
import json
from unittest import TestCase, main
from tempfile import NamedTemporaryFile

import numpy as np
import numpy.testing as npt

from biom.parse import (generatedby, MetadataMap, parse_biom_table, parse_uc,
                        load_table, save_table)
from biom.table import Table
from biom.util import __version__
from biom.tests.long_lines import (uc_empty, uc_invalid_id, uc_minimal,
                                   uc_lib_minimal,
                                   uc_seed_hits, uc_mixed_hits,
                                   uc_underscores_in_sample_id)
import h5py

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011-2017, The BIOM Format Development Team"
__credits__ = ["Justin Kuczynski", "Daniel McDonald", "Adam Robbins-Pianka",
               "Jose Antonio Navas Molina"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


class ParseTests(TestCase):

    """Tests of parse functions"""

    def setUp(self):
        """define some top-level data"""
        self.legacy_otu_table1 = legacy_otu_table1
        self.otu_table1 = otu_table1
        self.otu_table1_floats = otu_table1_floats
        self.to_remove = []
        self.biom_minimal_sparse = biom_minimal_sparse

        self.classic_otu_table1_w_tax = classic_otu_table1_w_tax.split('\n')
        self.classic_otu_table1_no_tax = classic_otu_table1_no_tax.split('\n')
        self.classic_table_with_complex_metadata = \
            classic_table_with_complex_metadata.split('\n')

    def tearDown(self):
        if self.to_remove:
            for f in self.to_remove:
                os.remove(f)

    def test_from_tsv_bug_854(self):
        data = StringIO('#FeatureID\tSample1')
        exp = Table([], [], ['Sample1'])
        obs = Table.from_tsv(data, None, None, lambda x: x)
        self.assertEqual(obs, exp)

    def test_generatedby(self):
        """get a generatedby string"""
        exp = "BIOM-Format %s" % __version__
        obs = generatedby()
        self.assertEqual(obs, exp)

    def test_metadata_map(self):
        """MetadataMap functions as expected

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        s1 = ['#sample\ta\tb', '#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        exp = {'x': {'a': 'y', 'b': 'z'}, 'i': {'a': 'j', 'b': 'k'}}
        obs = MetadataMap.from_file(s1)
        self.assertEqual(obs, exp)

        # check that we strip double quotes by default
        s2 = ['#sample\ta\tb', '#comment line to skip',
              '"x "\t" y "\t z ', ' ', '"#more skip"', 'i\t"j"\tk']
        obs = MetadataMap.from_file(s2)
        self.assertEqual(obs, exp)

    def test_metadata_map_w_map_fs(self):
        """MetadataMap functions as expected w process_fns

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        s1 = ['#sample\ta\tb', '#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        exp = {'x': {'a': 'y', 'b': 'zzz'}, 'i': {'a': 'j', 'b': 'kkk'}}
        process_fns = {'b': lambda x: x * 3}
        obs = MetadataMap.from_file(s1, process_fns=process_fns)
        self.assertEqual(obs, exp)

    def test_metadata_map_w_header(self):
        """MetadataMap functions as expected w user-provided header

        This method is ported from QIIME (http://www.qiime.org). QIIME is a GPL
        project, but we obtained permission from the authors of this method to
        port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
        """
        # number of user-provided headers matches number of columns, and no
        # header line in file
        s1 = ['#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        exp = {'x': {'a': 'y', 'b': 'z'}, 'i': {'a': 'j', 'b': 'k'}}
        header = ['sample', 'a', 'b']
        obs = MetadataMap.from_file(s1, header=header)
        self.assertEqual(obs, exp)

        # number of user-provided headers is fewer than number of columns, and
        # no header line in file
        s1 = ['#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a'],
               ['comment line to skip', 'more skip'])
        exp = {'x': {'a': 'y'}, 'i': {'a': 'j'}}
        header = ['sample', 'a']
        obs = MetadataMap.from_file(s1, header=header)
        self.assertEqual(obs, exp)

        # number of user-provided headers is fewer than number of columns, and
        # header line in file (overridden by user-provided)
        s1 = ['#sample\ta\tb', '#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a'],
               ['comment line to skip', 'more skip'])
        exp = {'x': {'a': 'y'}, 'i': {'a': 'j'}}
        header = ['sample', 'a']
        obs = MetadataMap.from_file(s1, header=header)
        self.assertEqual(obs, exp)

    def test_parse_biom_json(self):
        """test the biom otu table parser"""
        # light test. this code is used thoroughly within the other
        # parse_biom_table methods
        tab1_fh = json.load(StringIO(self.biom_minimal_sparse))
        tab = Table.from_json(tab1_fh)
        npt.assert_equal(tab.ids(), ('Sample1', 'Sample2', 'Sample3',
                                     'Sample4', 'Sample5', 'Sample6'))
        npt.assert_equal(tab.ids(axis='observation'),
                         ('GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3',
                          'GG_OTU_4', 'GG_OTU_5'))
        self.assertEqual(tab.metadata(), None)
        self.assertEqual(tab.metadata(axis='observation'), None)

        tab = parse_biom_table(StringIO(self.biom_minimal_sparse))
        npt.assert_equal(tab.ids(), ('Sample1', 'Sample2', 'Sample3',
                                     'Sample4', 'Sample5', 'Sample6'))
        npt.assert_equal(tab.ids(axis='observation'),
                         ('GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3',
                          'GG_OTU_4', 'GG_OTU_5'))
        self.assertEqual(tab.metadata(), None)
        self.assertEqual(tab.metadata(axis='observation'), None)

        tablestring = '''{
            "id":null,
            "format": "Biological Observation Matrix 0.9.1-dev",
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
        }'''
        tbs_fh = json.load(StringIO(tablestring))
        tab1 = Table.from_json(tbs_fh)
        tab2 = parse_biom_table(tablestring)
        self.assertEqual(tab1, tab2)

    def test_parse_biom_table_subset(self):
        """test the biom table parser subsetting"""
        tab = parse_biom_table(StringIO(self.biom_minimal_sparse),
                               ids=['Sample1', 'Sample3', 'Sample5',
                                    'Sample6'])
        npt.assert_equal(tab.ids(), ('Sample1', 'Sample3', 'Sample5',
                                     'Sample6'))
        npt.assert_equal(tab.ids(axis='observation'),
                         ('GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3', 'GG_OTU_4',
                          'GG_OTU_5'))
        self.assertEqual(tab.metadata(), None)
        self.assertEqual(tab.metadata(axis='observation'), None)

        tab = parse_biom_table(StringIO(self.biom_minimal_sparse),
                               ids=['GG_OTU_2', 'GG_OTU_3', 'GG_OTU_5'],
                               axis='observation')
        npt.assert_equal(tab.ids(), ('Sample1', 'Sample2', 'Sample3',
                                     'Sample4', 'Sample5', 'Sample6',))
        npt.assert_equal(tab.ids(axis='observation'),
                         ('GG_OTU_2', 'GG_OTU_3', 'GG_OTU_5'))
        self.assertEqual(tab.metadata(), None)
        self.assertEqual(tab.metadata(axis='observation'), None)

    def test_parse_adjacency_weird_input(self):
        with self.assertRaisesRegex(ValueError, "Not sure"):
            Table.from_adjacency({'foo', 'bar'})

    def test_parse_adjacency_bad_input(self):
        with self.assertRaisesRegex(ValueError, "Does not appear"):
            Table.from_adjacency(['a\tb\tc\n', 'd\te\tf\n'])

        with self.assertRaisesRegex(ValueError, "Does not appear"):
            Table.from_adjacency(['a\tb\n', 'd\te\n'])

        with self.assertRaises(AssertionError):
            Table.from_adjacency(['a\tb\t1\n', 'd\te\n'])

    def test_parse_adjacency_table_header(self):
        lines = ['#OTU ID\tSampleID\tvalue\n',
                 'O1\tS1\t10\n',
                 'O4\tS2\t1\n',
                 'O3\tS3\t2\n',
                 'O4\tS1\t5\n',
                 'O2\tS2\t3\n']
        exp = Table(np.array([[10, 0, 0],
                              [0, 3, 0],
                              [0, 0, 2],
                              [5, 1, 0]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3'])
        obs = Table.from_adjacency(''.join(lines))
        self.assertEqual(obs, exp)

    def test_parse_adjacency_table_no_header(self):
        lines = ['O1\tS1\t10\n',
                 'O4\tS2\t1\n',
                 'O3\tS3\t2\n',
                 'O4\tS1\t5\n',
                 'O2\tS2\t3\n']
        exp = Table(np.array([[10, 0, 0],
                              [0, 3, 0],
                              [0, 0, 2],
                              [5, 1, 0]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3'])
        obs = Table.from_adjacency(''.join(lines))
        self.assertEqual(obs, exp)

    def test_parse_biom_table_hdf5(self):
        """Make sure we can parse a HDF5 table through the same loader"""
        cwd = os.getcwd()
        if os.path.sep in __file__[1:]:
            os.chdir(os.path.dirname(__file__))
        Table.from_hdf5(h5py.File(os.path.join('test_data', 'test.biom'),
                                  'r'))
        os.chdir(cwd)

    def test_save_table_filepath(self):
        t = Table(np.array([[0, 1, 2], [3, 4, 5]]), ['a', 'b'],
                  ['c', 'd', 'e'])
        with NamedTemporaryFile(delete=False) as tmpfile:
            save_table(t, tmpfile.name)
            obs = load_table(tmpfile.name)
            self.assertEqual(obs, t)
        self.to_remove.append(tmpfile.name)

    def test_load_table_filepath(self):
        cwd = os.getcwd()
        if os.path.sep in __file__[1:]:
            os.chdir(os.path.dirname(__file__))
        load_table(os.path.join('test_data', 'test.biom'))
        os.chdir(cwd)

    def test_load_table_inmemory(self):
        cwd = os.getcwd()
        if os.path.sep in __file__[1:]:
            os.chdir(os.path.dirname(__file__))
        load_table(h5py.File(os.path.join('test_data', 'test.biom'), 'r'))
        os.chdir(cwd)

    def test_load_table_inmemory_json(self):
        cwd = os.getcwd()
        if os.path.sep in __file__[1:]:
            os.chdir(os.path.dirname(__file__))
        load_table(open(os.path.join('test_data', 'test.json')))
        os.chdir(cwd)

    def test_load_table_inmemory_stringio(self):
        load_table(StringIO('\n'.join(self.classic_otu_table1_no_tax)))

    def test_parse_biom_table(self):
        """tests for parse_biom_table when we do not have h5py"""
        # This is a TSV as a list of lines
        t = parse_biom_table(self.classic_otu_table1_no_tax)

        # Test TSV as a list of lines
        t_tsv_str = t.to_tsv()
        t_tsv_lines = t_tsv_str.splitlines()
        t_tsv = parse_biom_table(t_tsv_lines)
        self.assertEqual(t, t_tsv)

        # Test TSV as a file-like object
        t_tsv_stringio = StringIO(t_tsv_str)
        t_tsv = parse_biom_table(t_tsv_stringio)
        self.assertEqual(t, t_tsv)

        # Test JSON as a list of lines
        t_json_str = t.to_json('asd')
        t_json_lines = t_json_str.splitlines()
        t_json = parse_biom_table(t_json_lines)
        self.assertEqual(t, t_json)

        # Test JSON as a file-like object
        t_json_str = t.to_json('asd')
        t_json_stringio = StringIO(t_json_str)
        t_json = parse_biom_table(t_json_stringio)
        self.assertEqual(t, t_json)

    def test_parse_biom_table_with_hdf5(self):
        """tests for parse_biom_table when we have h5py"""
        # We will round-trip the HDF5 file to several different formats, and
        # make sure we can recover the same table using parse_biom_table
        if os.path.sep in __file__[1:]:
            os.chdir(os.path.dirname(__file__))

        t = parse_biom_table(h5py.File(os.path.join('test_data', 'test.biom'),
                                       'r'))

        # These things are not round-trippable using the general-purpose
        # parse_biom_table function
        t._sample_metadata = None
        t._observation_metadata = None
        t.type = None

        # Test TSV as a list of lines
        t_tsv_str = t.to_tsv()
        t_tsv_lines = t_tsv_str.splitlines()
        t_tsv = parse_biom_table(t_tsv_lines)
        self.assertEqual(t, t_tsv)
        # Test TSV as a file-like object
        t_tsv_stringio = StringIO(t_tsv_str)
        t_tsv = parse_biom_table(t_tsv_stringio)
        self.assertEqual(t, t_tsv)

        # Test JSON as a list of lines
        t_json_str = t.to_json('asd')
        t_json_lines = t_json_str.splitlines()
        t_json = parse_biom_table(t_json_lines)
        self.assertEqual(t, t_json)
        # Test JSON as a file-like object
        t_json_str = t.to_json('asd')
        t_json_stringio = StringIO(t_json_str)
        t_json = parse_biom_table(t_json_stringio)
        self.assertEqual(t, t_json)

    def test_empty_metadata_inconsistent_handling(self):
        oids = list('bacd')
        sids = list('YXZ')
        mat = np.array([[2, 1, 0], [0, 5, 0],
                        [0, 3, 0], [1, 2, 0]])

        A = Table(mat, oids, sids,
                  observation_metadata=[{}, {}, {}, {}],
                  sample_metadata=[{}, {}, {}])
        B = Table(mat, oids, sids)

        self.assertEqual(A, B)


legacy_otu_table1 = """# some comment goes here
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propioniba\
cterineae; Propionibacterium

1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; La\
ctobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniace\
ae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; St\
aphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1 = """# Some comment

OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propioniba\
cterineae; Propionibacterium
# some other comment
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; La\
ctobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniace\
ae; Corynebacteriaceae
# comments
#    everywhere!
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; St\
aphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1_floats = """# Some comment

OTU ID	Fing	Key	NA	Consensus Lineage
0	19111.0	44536.0	42.0	Bacteria; Actinobacteria; Actinobacteridae; Propio\
nibacterineae; Propionibacterium
# some other comment
1	1216.0	3500.0	6.0	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; La\
ctobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7	1803.0	1184.0	2.0	Bacteria; Actinobacteria; Actinobacteridae; Gordoniace\
ae; Corynebacteriaceae
# comments
#    everywhere!
3	1722.1	4903.2	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; St\
aphylococcaceae
4	589.6	2074.4	34.5	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

classic_table_with_complex_metadata = """# some comment
#OTU ID	sample1	sample2	KEGG_Pathways
K05842	1.0	3.5	rank1A; rank2A|rank1B; rank2B
K05841	2.0	4.5	Environmental Information Processing;
K00508	0.0	0.0	Metabolism; Lipid Metabolism; Linoleic acid metabolism
K00500	0.5	0.5	Metabolism; Amino Acid Metabolism; Phenylalanine metabolism|Me\
tabolism; Amino Acid Metabolism; Phenylalanine, tyrosine and tryptophan biosyn\
thesis
K00507	0.0	0.0	Metabolism; Lipid Metabolism; Biosynthesis of unsaturated fatt\
y acids|Organismal Systems; Endocrine System; PPAR signaling pathway
"""

biom_minimal_sparse = """
    {
        "id":null,
        "format": "Biological Observation Matrix v0.9",
        "format_url": "http://some_website/QIIME_MGRAST_dataformat_v0.9.html",
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
                [2,4,2],
                [3,0,2],
                [3,1,1],
                [3,2,1],
                [3,5,1],
                [4,1,1],
                [4,2,1]
               ]
    }
"""

classic_otu_table1_w_tax = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636\
\tConsensus Lineage
0	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
1	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
2	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Porphyromonadaceae;Parabacteroides
3	2	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
4	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
5	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
6	0	0	0	0	0	0	0	1	0	Root;Bacteria;Actinobacteria;Actinobac\
teria
7	0	0	2	0	0	0	0	0	2	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
8	1	1	0	2	4	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lacto\
bacillales;Lactobacillaceae;Lactobacillus
9	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
10	0	1	0	0	0	0	0	0	0	Root;Bacteria
11	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
12	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
13	1	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
14	0	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
15	0	0	0	0	1	0	0	0	0	Root;Bacteria
16	1	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
17	0	0	0	1	0	0	4	10	37	Root;Bacteria;Bacteroidetes
18	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
19	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
20	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
21	0	0	0	0	0	0	2	3	2	Root;Bacteria;Bacteroidetes
22	0	0	0	0	2	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
23	14	1	14	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lacto\
bacillales;Lactobacillaceae;Lactobacillus
24	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
25	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
26	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
27	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
28	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
29	6	0	4	0	2	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
30	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
31	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
32	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
33	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
34	0	0	0	0	0	0	8	10	2	Root;Bacteria
35	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
36	1	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
37	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
38	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
39	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroid\
etes;Bacteroidales;Rikenellaceae;Alistipes
40	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
41	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
42	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
43	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
44	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
45	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Erysipelotric\
hi;Erysipelotrichales;Erysipelotrichaceae;Coprobacillus
46	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
47	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
48	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
49	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
50	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
51	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
52	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
53	0	0	0	0	0	0	2	0	1	Root;Bacteria;Proteobacteria;Deltaprot\
eobacteria
54	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Porphyromonadaceae;Parabacteroides
55	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
56	0	0	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes
57	0	0	0	0	0	0	0	1	0	Root;Bacteria
58	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
59	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferrib\
acteres;Deferribacterales;Deferribacteraceae;Mucispirillum
60	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
61	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
62	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
63	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
64	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
65	0	0	0	6	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
66	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
67	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
68	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
69	0	0	1	0	0	0	0	0	0	Root;Bacteria
70	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
71	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
72	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
73	0	0	0	0	0	5	0	0	0	Root;Bacteria;Bacteroidetes
74	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
75	1	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
76	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
77	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
78	1	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
79	2	3	8	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
80	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Porphyromonadaceae;Parabacteroides
81	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
82	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
83	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
84	1	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae;Ruminococcus
85	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
86	0	0	0	0	0	0	0	1	0	Root;Bacteria
87	0	0	1	0	0	2	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
88	0	0	0	0	0	0	0	1	0	Root;Bacteria
89	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
90	0	0	0	9	0	0	3	0	0	Root;Bacteria;Firmicutes;Erysipelotric\
hi;Erysipelotrichales;Erysipelotrichaceae;Turicibacter
91	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Butyrivibrio
92	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
93	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
94	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
95	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
96	0	0	0	1	0	1	0	1	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
97	0	0	0	0	0	1	0	0	0	Root;Bacteria
98	0	0	0	0	0	0	0	1	0	Root;Bacteria
99	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
100	0	0	0	1	0	0	0	0	0	Root;Bacteria
101	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
102	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
103	0	1	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
104	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
105	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
106	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
107	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
108	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Incertae Sedis XIII;Anaerovorax
109	0	0	0	1	0	0	1	5	2	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
110	0	0	0	0	0	2	0	0	0	Root;Bacteria;Actinobacteria;Actinobac\
teria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Olse\
nella
111	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
112	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
113	0	0	0	0	0	1	0	0	0	Root;Bacteria
114	0	0	0	0	0	1	0	0	0	Root;Bacteria
115	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
116	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
117	1	0	2	0	0	6	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
118	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
119	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
120	1	3	1	2	1	9	2	4	5	Root;Bacteria;Bacteroidetes
121	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
122	0	0	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
123	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobac\
teria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
124	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobac\
teria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
125	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
126	0	0	2	0	0	0	0	1	0	Root;Bacteria
127	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
128	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
129	0	0	0	1	0	0	0	0	0	Root;Bacteria
130	0	0	0	0	5	2	0	0	0	Root;Bacteria;Proteobacteria;Epsilonpr\
oteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter
131	0	0	1	3	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
132	0	0	0	0	1	0	0	0	0	Root;Bacteria
133	0	0	1	0	0	0	0	0	0	Root;Bacteria
134	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
135	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
136	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
137	0	0	0	0	0	0	0	1	0	Root;Bacteria
138	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
139	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
140	0	0	0	0	0	0	1	3	0	Root;Bacteria
141	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
142	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
143	0	0	1	0	0	0	0	0	0	Root;Bacteria
144	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
145	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
146	1	0	0	0	2	0	2	0	3	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
147	0	1	0	1	1	0	0	0	3	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
148	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
149	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
150	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
151	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
152	0	0	0	1	0	0	1	2	19	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
153	0	2	1	2	0	0	1	1	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
154	2	18	0	1	0	0	21	4	4	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
155	0	0	0	0	0	5	9	5	3	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
156	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
157	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
158	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
159	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
160	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
161	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
162	0	0	0	0	0	3	5	2	6	Root;Bacteria;Deferribacteres;Deferrib\
acteres;Deferribacterales;Deferribacteraceae;Mucispirillum
163	0	0	0	0	0	0	0	0	1	Root;Bacteria
164	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
165	2	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
166	0	0	0	0	0	0	0	1	0	Root;Bacteria
167	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
168	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
169	0	2	0	7	0	0	0	2	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
170	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
171	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
172	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
173	0	0	0	0	0	1	0	0	0	Root;Bacteria
174	1	0	0	0	10	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Peptostreptococcaceae;Peptostreptococcaceae Incertae Sedis
175	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
176	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
177	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia
178	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
179	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
180	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
181	1	4	2	6	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
182	0	0	0	0	0	1	0	0	0	Root;Bacteria
183	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia
184	0	0	0	1	0	0	3	1	0	Root;Bacteria;Bacteroidetes
185	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
186	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
187	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
188	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
189	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
190	0	0	0	0	0	0	0	1	0	Root;Bacteria
191	2	1	10	2	24	0	0	1	1	Root;Bacteria
192	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lacto\
bacillales;Streptococcaceae;Streptococcus
193	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Butyrivibrio
194	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae;Acetanaerobacterium
195	0	0	0	0	0	1	0	0	0	Root;Bacteria
196	0	0	0	0	0	1	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
197	0	1	0	0	0	0	0	0	0	Root;Bacteria
198	0	2	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales
199	0	0	0	0	0	1	1	0	0	Root;Bacteria
200	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
201	0	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
202	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
203	0	2	2	4	0	5	1	5	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
204	1	4	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
205	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
206	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
207	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
208	0	2	0	2	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
209	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
210	0	0	0	0	0	0	0	0	1	Root;Bacteria
211	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
212	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
213	0	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes
214	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
215	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
216	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
217	0	0	0	0	0	2	0	1	0	Root;Bacteria
218	0	0	0	0	9	1	0	0	0	Root;Bacteria;Bacteroidetes
219	0	0	0	0	1	0	0	0	0	Root;Bacteria
220	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
221	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
222	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
223	0	0	0	0	0	0	0	2	2	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
224	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
225	0	2	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
226	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
227	0	1	2	0	9	1	1	1	3	Root;Bacteria;Bacteroidetes
228	16	0	0	0	12	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
229	0	0	0	0	0	1	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Incertae Sedis XIII
230	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
231	0	19	2	0	2	0	3	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
232	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
233	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
234	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lacto\
bacillales;Lactobacillaceae;Lactobacillus
235	0	1	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
236	0	0	0	0	0	2	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales
237	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
238	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
239	0	0	0	0	0	1	0	0	0	Root;Bacteria
240	0	0	0	0	0	1	0	0	0	Root;Bacteria
241	0	0	0	0	0	0	2	0	0	Root;Bacteria;TM7;TM7_genera_incertae_\
sedis
242	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
243	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
244	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
245	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
246	0	0	0	0	0	0	0	1	0	Root;Bacteria
247	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
248	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lacto\
bacillales;Lactobacillaceae;Lactobacillus
249	1	0	0	0	0	0	0	0	0	Root;Bacteria
250	1	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
251	0	0	0	1	4	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
252	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
253	0	0	0	0	2	0	0	5	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
254	11	13	6	13	2	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
255	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
256	0	0	0	0	0	0	1	0	0	Root;Bacteria
257	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
258	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
259	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
260	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
261	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
262	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Bryantella
263	0	0	0	0	1	0	0	0	0	Root;Bacteria
264	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
265	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
266	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae;Alistipes
267	1	0	0	5	17	20	0	0	0	Root;Bacteria
268	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
269	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
270	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
271	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
272	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
273	0	0	0	0	0	0	1	0	0	Root;Bacteria
274	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
275	0	0	0	0	0	0	1	0	0	Root;Bacteria;Verrucomicrobia;Verrucom\
icrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
276	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
277	1	0	0	0	0	0	0	0	0	Root;Bacteria
278	0	0	0	0	0	1	0	0	0	Root;Bacteria
279	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
280	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
281	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
282	0	0	0	0	0	0	2	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Porphyromonadaceae;Parabacteroides
283	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
284	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
285	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
286	0	2	3	1	4	0	5	0	4	Root;Bacteria;Bacteroidetes
287	0	0	0	0	0	0	1	1	1	Root;Bacteria;Bacteroidetes
288	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
289	0	0	0	0	3	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
290	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;Bacilli;Bacil\
lales;Staphylococcaceae;Staphylococcus
291	0	0	0	0	1	0	0	0	0	Root;Bacteria
292	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
293	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
294	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
295	29	1	10	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
296	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
297	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
298	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobac\
teria
299	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
300	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia
301	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
302	0	0	0	0	0	1	0	0	0	Root;Bacteria
303	0	0	0	0	0	0	0	0	1	Root;Bacteria
304	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
305	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
306	0	0	0	0	0	0	0	0	1	Root;Bacteria
307	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
308	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae;Ruminococcaceae Incertae Sedis
309	0	0	0	1	0	0	0	0	0	Root;Bacteria;Actinobacteria;Actinobac\
teria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Deni\
trobacterium
310	0	0	1	0	0	0	0	0	0	Root;Bacteria
311	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
312	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
313	0	1	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Porphyromonadaceae;Parabacteroides
314	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
315	1	3	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
316	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
317	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
318	0	0	0	0	0	1	0	0	0	Root;Bacteria;Proteobacteria
319	0	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
320	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
321	0	0	0	0	0	0	0	0	1	Root;Bacteria
322	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
323	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
324	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
325	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
326	0	0	0	0	4	0	0	0	2	Root;Bacteria;Firmicutes;Erysipelotric\
hi;Erysipelotrichales;Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
327	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
328	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
329	2	2	0	1	0	0	0	0	0	Root;Bacteria;Bacteroidetes
330	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes
331	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes
332	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
333	0	0	0	0	0	6	0	3	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
334	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
335	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
336	0	0	1	0	0	0	0	0	0	Root;Bacteria
337	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
338	0	0	0	0	0	0	0	1	0	Root;Bacteria
339	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
340	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
341	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
342	0	0	0	0	0	1	0	0	0	Root;Bacteria
343	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobac\
teria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
344	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
345	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
346	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
347	0	0	0	1	0	0	0	0	0	Root;Bacteria
348	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
349	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
350	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
351	0	0	0	0	2	2	1	4	1	Root;Bacteria;Bacteroidetes
352	3	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
353	0	4	4	0	1	2	0	2	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
354	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
355	0	0	0	0	0	0	0	1	0	Root;Bacteria
356	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
357	0	0	0	4	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
358	0	0	1	0	0	0	0	0	0	Root;Bacteria
359	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
360	0	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
361	2	0	2	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
362	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
363	0	0	0	0	0	1	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Rikenellaceae
364	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
365	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
366	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Roseburia
367	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
368	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
369	0	0	0	0	0	1	0	0	0	Root;Bacteria
370	2	1	0	5	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
371	1	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
372	0	1	0	0	0	0	0	0	0	Root;Bacteria
373	0	1	0	0	0	0	3	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Clostridiaceae;Clostridiaceae 1;Clostridium
374	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
375	0	0	0	0	0	0	4	0	0	Root;Bacteria;Firmicutes;Erysipelotric\
hi;Erysipelotrichales;Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
376	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
377	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
378	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
379	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Ruminococcaceae
380	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Bacilli;Bacil\
lales;Staphylococcaceae;Staphylococcus
381	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
382	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
383	4	9	0	2	0	0	0	2	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
384	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
385	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Bacilli;Lacto\
bacillales;Carnobacteriaceae;Carnobacteriaceae 1
386	0	0	1	0	0	0	0	0	0	Root;Bacteria
387	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
388	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
389	0	1	0	0	0	0	0	0	0	Root;Bacteria
390	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
391	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes
392	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
393	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
394	0	0	1	0	0	0	0	0	0	Root;Bacteria
395	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
396	2	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
397	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
398	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
399	0	0	0	0	0	0	13	0	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Bacteroidaceae;Bacteroides
400	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
401	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
402	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
403	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroide\
tes;Bacteroidales;Prevotellaceae
404	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
405	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
406	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
407	1	0	0	0	0	4	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
408	1	5	3	2	0	0	0	0	1	Root;Bacteria;Bacteroidetes
409	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
410	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
411	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
412	0	0	0	0	2	0	0	0	0	Root;Bacteria;Bacteroidetes
413	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales
414	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales;Lachnospiraceae
415	0	0	0	0	0	7	0	2	2	Root;Bacteria;Bacteroidetes
416	0	1	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Cl\
ostridiales"""

classic_otu_table1_no_tax = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
0	0	0	0	0	0	0	0	1	0
1	0	0	0	0	0	1	0	0	0
2	0	0	0	0	0	0	0	0	1
3	2	1	0	0	0	0	0	0	0
4	1	0	0	0	0	0	0	0	0
5	0	0	0	0	0	0	0	0	1
6	0	0	0	0	0	0	0	1	0
7	0	0	2	0	0	0	0	0	2
8	1	1	0	2	4	0	0	0	0
9	0	0	2	0	0	0	0	0	0
10	0	1	0	0	0	0	0	0	0
11	0	0	0	0	0	0	1	0	0
12	0	0	0	0	0	0	1	0	0
13	1	0	0	1	0	1	0	0	0
14	0	0	1	1	0	0	0	0	0
15	0	0	0	0	1	0	0	0	0
16	1	0	2	0	0	0	0	0	0
17	0	0	0	1	0	0	4	10	37
18	0	1	0	0	0	0	0	0	0
19	0	0	0	0	0	0	0	0	1
20	0	0	0	0	1	0	0	0	0
21	0	0	0	0	0	0	2	3	2
22	0	0	0	0	2	0	1	0	0
23	14	1	14	1	0	0	0	0	0
24	1	0	0	0	0	0	0	0	0
25	0	0	0	1	0	0	0	0	0
26	0	0	0	0	0	0	0	1	1
27	0	0	0	0	0	0	0	0	1
28	0	1	0	0	0	0	0	0	0
29	6	0	4	0	2	0	0	0	0
30	0	0	0	0	0	1	0	0	0
31	1	0	0	0	0	0	0	0	0
32	0	0	0	0	1	0	0	0	0
33	0	0	0	1	0	0	0	0	0
34	0	0	0	0	0	0	8	10	2
35	1	0	1	0	0	0	0	0	0
36	1	0	1	0	0	0	0	1	1
37	0	0	0	0	0	1	0	0	0
38	0	0	1	0	0	0	0	0	0
39	0	0	0	0	0	0	0	1	0
40	0	0	1	0	0	0	0	0	0
41	0	0	1	0	0	0	0	1	0
42	0	0	0	0	0	1	0	0	0
43	0	0	0	0	0	1	0	0	0
44	0	0	1	0	0	0	0	0	0
45	1	0	0	0	0	0	0	0	0
46	0	0	0	0	0	0	0	0	1
47	0	0	0	1	0	0	0	0	0
48	0	0	0	0	1	0	0	0	0
49	0	0	0	1	0	0	0	0	0
50	0	1	0	0	0	0	0	0	0
51	0	1	0	0	0	0	0	0	0
52	0	2	0	0	0	0	0	0	0
53	0	0	0	0	0	0	2	0	1
54	0	0	0	0	0	0	5	0	0
55	0	0	0	0	0	0	1	0	0
56	0	0	0	0	0	1	0	0	0
57	0	0	0	0	0	0	0	1	0
58	1	0	1	0	0	0	0	0	0
59	0	0	0	0	0	0	0	0	1
60	0	0	0	0	0	0	0	1	0
61	0	0	1	0	0	0	0	1	0
62	0	0	1	0	0	0	0	0	0
63	1	0	1	0	0	0	0	0	0
64	0	0	0	0	0	0	0	0	1
65	0	0	0	6	0	0	0	1	0
66	0	0	1	0	0	0	0	0	0
67	0	0	1	0	0	0	0	0	0
68	1	0	0	0	0	0	0	0	0
69	0	0	1	0	0	0	0	0	0
70	0	0	0	0	0	1	0	0	0
71	0	0	1	0	0	0	0	0	0
72	0	0	0	0	0	1	0	0	0
73	0	0	0	0	0	5	0	0	0
74	0	0	0	1	0	0	0	0	0
75	1	0	1	0	0	0	0	0	0
76	0	0	0	1	0	0	0	0	0
77	0	0	0	1	0	0	0	0	0
78	1	0	1	1	0	0	0	0	0
79	2	3	8	0	1	0	0	0	0
80	0	0	0	0	0	0	0	0	1
81	1	0	0	0	0	0	0	0	0
82	0	0	0	0	0	2	0	0	0
83	0	0	0	1	0	0	0	1	0
84	1	0	0	0	0	0	0	2	0
85	0	0	0	0	0	0	0	0	1
86	0	0	0	0	0	0	0	1	0
87	0	0	1	0	0	2	0	1	0
88	0	0	0	0	0	0	0	1	0
89	0	0	1	0	0	0	0	0	0
90	0	0	0	9	0	0	3	0	0
91	0	0	0	1	0	0	0	0	0
92	0	0	0	0	0	0	1	0	0
93	0	0	0	0	0	0	2	1	0
94	0	0	0	0	0	0	0	1	0
95	0	0	0	2	0	0	0	0	0
96	0	0	0	1	0	1	0	1	1
97	0	0	0	0	0	1	0	0	0
98	0	0	0	0	0	0	0	1	0
99	0	0	0	1	0	0	0	0	0
100	0	0	0	1	0	0	0	0	0
101	0	0	0	3	0	0	0	0	0
102	0	1	0	0	0	0	0	0	0
103	0	1	0	0	0	0	1	0	0
104	0	0	0	0	0	1	0	0	0
105	0	1	0	0	0	0	0	0	0
106	0	0	0	0	0	1	0	0	0
107	0	0	0	0	0	1	0	0	0
108	0	0	0	0	0	0	1	0	0
109	0	0	0	1	0	0	1	5	2
110	0	0	0	0	0	2	0	0	0
111	0	0	0	0	0	0	1	0	0
112	0	0	0	0	0	0	1	0	0
113	0	0	0	0	0	1	0	0	0
114	0	0	0	0	0	1	0	0	0
115	0	0	0	0	0	1	0	0	0
116	0	1	0	0	0	0	0	0	0
117	1	0	2	0	0	6	0	0	0
118	0	0	0	1	0	0	0	0	0
119	0	0	0	0	0	0	0	1	0
120	1	3	1	2	1	9	2	4	5
121	0	0	0	0	0	0	0	1	0
122	0	0	0	1	0	2	0	0	0
123	0	0	0	0	0	0	1	0	0
124	0	0	0	0	0	0	1	0	0
125	0	0	0	0	0	0	1	0	0
126	0	0	2	0	0	0	0	1	0
127	0	0	0	0	0	1	0	0	0
128	0	0	0	0	0	0	1	0	0
129	0	0	0	1	0	0	0	0	0
130	0	0	0	0	5	2	0	0	0
131	0	0	1	3	0	0	0	0	0
132	0	0	0	0	1	0	0	0	0
133	0	0	1	0	0	0	0	0	0
134	0	0	0	0	0	0	0	0	1
135	0	0	1	0	0	0	0	0	0
136	1	0	0	0	0	0	0	0	0
137	0	0	0	0	0	0	0	1	0
138	0	0	1	0	0	0	0	0	0
139	1	0	0	0	0	0	0	0	0
140	0	0	0	0	0	0	1	3	0
141	0	0	0	0	1	0	0	0	0
142	0	0	0	0	1	0	0	0	0
143	0	0	1	0	0	0	0	0	0
144	0	0	0	0	0	1	0	0	0
145	0	0	2	0	0	0	0	0	0
146	1	0	0	0	2	0	2	0	3
147	0	1	0	1	1	0	0	0	3
148	0	0	0	0	0	1	0	0	0
149	0	0	0	0	0	0	0	1	0
150	0	0	0	0	1	0	0	0	0
151	0	0	0	1	0	0	0	1	0
152	0	0	0	1	0	0	1	2	19
153	0	2	1	2	0	0	1	1	1
154	2	18	0	1	0	0	21	4	4
155	0	0	0	0	0	5	9	5	3
156	0	0	1	0	0	0	0	1	0
157	0	0	1	0	0	0	0	0	0
158	1	0	1	0	0	0	0	0	0
159	0	0	0	0	0	0	0	1	1
160	0	0	0	0	0	0	1	0	0
161	0	0	1	0	0	0	0	0	0
162	0	0	0	0	0	3	5	2	6
163	0	0	0	0	0	0	0	0	1
164	0	0	0	0	0	1	0	0	0
165	2	1	1	0	0	0	0	0	0
166	0	0	0	0	0	0	0	1	0
167	1	0	0	0	0	0	0	0	0
168	0	0	0	1	0	0	0	0	0
169	0	2	0	7	0	0	0	2	0
170	0	0	0	1	0	0	0	0	0
171	0	0	0	1	0	0	0	0	0
172	1	0	0	0	0	0	0	0	0
173	0	0	0	0	0	1	0	0	0
174	1	0	0	0	10	0	0	0	0
175	0	0	0	0	1	0	0	0	0
176	0	0	0	0	0	1	0	0	0
177	0	0	0	1	0	0	0	0	0
178	0	0	0	2	0	0	0	0	0
179	0	0	0	1	0	0	0	0	0
180	0	0	0	0	1	0	0	0	0
181	1	4	2	6	0	0	0	0	0
182	0	0	0	0	0	1	0	0	0
183	0	0	0	0	0	0	1	0	0
184	0	0	0	1	0	0	3	1	0
185	0	0	0	0	0	0	0	0	1
186	0	0	1	0	0	0	0	0	0
187	0	1	0	0	0	0	0	0	1
188	0	0	0	0	0	0	0	1	0
189	0	0	0	1	0	0	0	0	0
190	0	0	0	0	0	0	0	1	0
191	2	1	10	2	24	0	0	1	1
192	0	0	0	0	0	1	0	0	0
193	0	0	0	0	0	1	0	0	0
194	0	0	2	0	0	0	0	0	0
195	0	0	0	0	0	1	0	0	0
196	0	0	0	0	0	1	0	1	0
197	0	1	0	0	0	0	0	0	0
198	0	2	0	0	0	1	0	0	0
199	0	0	0	0	0	1	1	0	0
200	0	0	0	2	0	0	0	0	0
201	0	0	0	1	0	1	0	0	0
202	0	0	0	0	0	0	1	0	0
203	0	2	2	4	0	5	1	5	0
204	1	4	0	1	0	0	0	0	0
205	0	0	0	0	0	0	0	1	0
206	0	1	0	0	0	0	0	0	0
207	0	0	0	0	0	0	0	1	0
208	0	2	0	2	0	0	0	1	0
209	0	0	1	0	0	0	0	0	0
210	0	0	0	0	0	0	0	0	1
211	1	0	0	1	0	0	0	0	0
212	0	0	0	0	0	0	0	0	1
213	0	0	0	0	0	0	0	2	0
214	0	0	0	0	0	0	0	1	0
215	0	0	0	0	0	0	0	1	0
216	0	0	0	0	0	0	0	1	0
217	0	0	0	0	0	2	0	1	0
218	0	0	0	0	9	1	0	0	0
219	0	0	0	0	1	0	0	0	0
220	1	0	0	0	1	0	0	0	0
221	0	0	0	0	0	0	0	1	0
222	0	1	0	0	0	0	0	0	0
223	0	0	0	0	0	0	0	2	2
224	0	0	0	1	0	0	0	0	0
225	0	2	1	0	0	0	0	0	0
226	0	0	0	0	0	1	0	0	0
227	0	1	2	0	9	1	1	1	3
228	16	0	0	0	12	0	0	0	0
229	0	0	0	0	0	1	1	0	0
230	0	0	0	1	0	0	0	0	0
231	0	19	2	0	2	0	3	0	0
232	0	0	0	0	0	0	1	0	0
233	0	0	0	0	1	0	0	0	0
234	0	0	0	0	1	0	0	0	0
235	0	1	1	0	1	0	0	0	0
236	0	0	0	0	0	2	0	0	0
237	0	0	0	0	1	0	0	0	0
238	0	0	0	0	0	0	0	1	0
239	0	0	0	0	0	1	0	0	0
240	0	0	0	0	0	1	0	0	0
241	0	0	0	0	0	0	2	0	0
242	0	0	0	0	0	0	1	0	0
243	0	0	0	0	0	0	1	0	0
244	0	0	0	0	0	0	0	0	1
245	0	0	0	1	0	0	0	1	0
246	0	0	0	0	0	0	0	1	0
247	0	0	1	0	0	0	0	0	0
248	1	0	0	1	0	0	0	0	0
249	1	0	0	0	0	0	0	0	0
250	1	0	0	0	0	0	0	1	0
251	0	0	0	1	4	0	0	0	0
252	0	0	0	1	0	0	0	0	0
253	0	0	0	0	2	0	0	5	0
254	11	13	6	13	2	0	0	0	0
255	0	0	0	0	0	1	0	0	0
256	0	0	0	0	0	0	1	0	0
257	0	0	0	0	0	0	5	0	0
258	0	0	1	0	0	0	0	0	0
259	0	0	0	0	0	0	0	1	0
260	0	0	0	0	0	0	0	1	0
261	0	0	0	0	0	0	0	1	0
262	0	1	0	0	0	0	0	0	1
263	0	0	0	0	1	0	0	0	0
264	0	0	0	0	0	1	0	0	0
265	0	0	0	0	0	2	0	0	0
266	0	0	0	2	0	0	0	0	0
267	1	0	0	5	17	20	0	0	0
268	0	0	0	0	0	0	1	0	0
269	0	0	0	1	0	0	0	0	0
270	0	0	1	0	0	0	0	0	0
271	0	0	0	0	0	0	0	0	1
272	0	0	0	1	0	0	0	0	0
273	0	0	0	0	0	0	1	0	0
274	0	0	0	0	0	0	1	0	0
275	0	0	0	0	0	0	1	0	0
276	0	0	0	0	0	0	0	1	0
277	1	0	0	0	0	0	0	0	0
278	0	0	0	0	0	1	0	0	0
279	0	0	0	0	0	1	0	0	0
280	0	1	0	0	0	0	0	0	0
281	1	0	0	0	0	0	0	0	0
282	0	0	0	0	0	0	2	0	0
283	0	0	0	0	0	0	2	1	0
284	0	0	0	1	0	0	0	0	0
285	0	0	0	0	0	0	1	0	0
286	0	2	3	1	4	0	5	0	4
287	0	0	0	0	0	0	1	1	1
288	0	0	0	0	0	1	0	0	0
289	0	0	0	0	3	0	0	0	0
290	0	0	0	0	0	0	0	0	2
291	0	0	0	0	1	0	0	0	0
292	0	0	0	0	1	0	0	0	0
293	0	0	0	0	0	1	0	0	0
294	0	1	0	0	0	0	0	0	0
295	29	1	10	0	0	0	0	0	0
296	0	0	0	0	1	0	0	0	0
297	0	0	0	1	0	0	0	0	0
298	0	0	0	0	0	0	1	0	0
299	0	0	0	0	0	0	1	0	1
300	0	0	0	0	0	1	0	0	0
301	0	0	0	0	0	0	2	0	0
302	0	0	0	0	0	1	0	0	0
303	0	0	0	0	0	0	0	0	1
304	0	0	0	0	0	0	0	1	0
305	1	0	0	0	0	0	0	0	0
306	0	0	0	0	0	0	0	0	1
307	0	0	1	0	0	0	0	0	0
308	0	1	0	0	0	0	0	0	0
309	0	0	0	1	0	0	0	0	0
310	0	0	1	0	0	0	0	0	0
311	0	0	0	0	0	1	0	0	0
312	0	0	1	0	0	0	0	0	0
313	0	1	0	0	0	0	0	0	1
314	0	0	1	0	0	0	0	0	0
315	1	3	1	0	0	0	0	0	0
316	0	1	0	0	0	0	0	0	0
317	0	0	0	0	0	0	1	0	0
318	0	0	0	0	0	1	0	0	0
319	0	2	1	0	0	0	0	0	0
320	0	0	0	1	0	0	0	0	0
321	0	0	0	0	0	0	0	0	1
322	0	0	0	1	0	0	0	0	0
323	0	0	1	0	0	0	0	0	0
324	0	0	1	0	0	0	0	0	0
325	0	1	0	0	0	0	0	0	0
326	0	0	0	0	4	0	0	0	2
327	0	0	0	0	0	0	0	1	0
328	0	0	0	1	0	0	0	0	0
329	2	2	0	1	0	0	0	0	0
330	0	0	1	0	0	0	0	0	0
331	0	0	0	0	1	0	0	0	0
332	0	1	0	0	0	0	0	0	0
333	0	0	0	0	0	6	0	3	0
334	1	0	0	0	0	0	0	0	0
335	0	0	0	0	0	0	0	1	0
336	0	0	1	0	0	0	0	0	0
337	0	0	0	1	0	0	0	0	0
338	0	0	0	0	0	0	0	1	0
339	0	0	1	0	0	0	0	0	0
340	0	0	2	0	0	0	0	0	0
341	0	0	1	0	0	0	0	1	0
342	0	0	0	0	0	1	0	0	0
343	0	0	0	0	0	0	0	0	1
344	0	0	1	0	0	0	0	0	0
345	1	0	0	0	0	0	0	0	0
346	0	1	0	0	0	0	0	0	0
347	0	0	0	1	0	0	0	0	0
348	0	0	0	1	0	0	0	0	0
349	0	0	0	0	0	0	1	0	1
350	1	0	0	0	0	0	0	0	0
351	0	0	0	0	2	2	1	4	1
352	3	0	0	0	0	0	0	0	0
353	0	4	4	0	1	2	0	2	1
354	0	0	0	0	0	1	0	0	0
355	0	0	0	0	0	0	0	1	0
356	0	0	0	0	0	1	0	0	0
357	0	0	0	4	0	0	0	0	0
358	0	0	1	0	0	0	0	0	0
359	0	0	1	0	0	0	0	0	0
360	0	0	1	0	0	0	0	1	1
361	2	0	2	1	0	0	0	0	0
362	1	0	0	1	0	0	0	0	0
363	0	0	0	0	0	1	0	1	0
364	1	0	0	0	0	0	0	0	0
365	0	0	0	0	0	2	0	0	0
366	0	0	0	1	0	0	0	0	0
367	0	0	0	0	1	0	0	0	0
368	0	0	0	0	0	1	0	0	0
369	0	0	0	0	0	1	0	0	0
370	2	1	0	5	0	1	0	0	0
371	1	1	0	0	0	0	0	0	0
372	0	1	0	0	0	0	0	0	0
373	0	1	0	0	0	0	3	0	0
374	0	0	0	0	0	0	1	0	0
375	0	0	0	0	0	0	4	0	0
376	0	0	0	0	0	0	0	0	1
377	0	0	0	0	0	0	0	1	0
378	0	0	0	0	0	0	0	0	1
379	0	0	0	0	0	1	0	0	0
380	0	0	0	0	0	0	0	0	1
381	0	0	2	0	0	0	0	0	0
382	0	0	0	0	0	0	0	1	0
383	4	9	0	2	0	0	0	2	0
384	0	0	1	0	0	0	0	0	0
385	0	0	0	0	0	0	0	0	1
386	0	0	1	0	0	0	0	0	0
387	0	0	1	0	0	0	0	0	0
388	0	0	0	1	0	0	0	0	0
389	0	1	0	0	0	0	0	0	0
390	0	0	0	0	0	0	0	0	1
391	0	0	0	0	0	0	0	0	1
392	0	1	0	0	0	0	0	0	0
393	0	0	0	0	0	1	0	0	0
394	0	0	1	0	0	0	0	0	0
395	1	1	1	0	0	0	0	0	0
396	2	0	0	0	0	0	0	0	0
397	0	0	0	0	0	0	0	1	0
398	0	0	0	0	0	0	0	1	0
399	0	0	0	0	0	0	13	0	0
400	0	0	0	0	0	0	1	0	0
401	0	1	0	0	0	0	0	0	0
402	0	1	0	0	0	0	0	0	0
403	0	0	0	0	0	0	0	1	0
404	0	0	0	0	0	0	0	1	0
405	0	0	0	0	0	0	0	1	0
406	0	0	0	0	0	1	0	0	0
407	1	0	0	0	0	4	0	0	0
408	1	5	3	2	0	0	0	0	1
409	0	0	0	0	0	0	0	1	1
410	0	0	0	0	1	0	0	0	0
411	0	0	0	1	0	0	0	0	0
412	0	0	0	0	2	0	0	0	0
413	0	0	0	0	0	0	0	1	0
414	1	0	1	0	0	0	0	0	0
415	0	0	0	0	0	7	0	2	2
416	0	1	0	0	1	0	0	0	0"""


class ParseUcTests(TestCase):

    def test_empty(self):
        """ empty uc file returns empty Table
        """
        actual = parse_uc(uc_empty.split('\n'))
        expected = Table(np.array([[]]),
                         observation_ids=[],
                         sample_ids=[])
        self.assertEqual(actual, expected)

    def test_minimal(self):
        """ single new seed observed
        """
        actual = parse_uc(uc_minimal.split('\n'))
        expected = Table(np.array([[1.0]]),
                         observation_ids=['f2_1539'],
                         sample_ids=['f2'])
        self.assertEqual(actual, expected)

    def test_lib_minimal(self):
        """ single library seed observed
        """
        actual = parse_uc(uc_lib_minimal.split('\n'))
        expected = Table(np.array([[1.0]]),
                         observation_ids=['295053'],
                         sample_ids=['f2'])
        self.assertEqual(actual, expected)

    def test_invalid(self):
        """ invalid query sequence identifier detected
        """
        self.assertRaises(ValueError, parse_uc, uc_invalid_id.split('\n'))

    def test_seed_hits(self):
        """ multiple new seeds observed
        """
        actual = parse_uc(uc_seed_hits.split('\n'))
        expected = Table(np.array([[2.0, 1.0], [0.0, 1.0]]),
                         observation_ids=['f2_1539', 'f3_44'],
                         sample_ids=['f2', 'f3'])
        self.assertEqual(actual, expected)

    def test_mixed_hits(self):
        """ new and library seeds observed
        """
        actual = parse_uc(uc_mixed_hits.split('\n'))
        expected = Table(np.array([[2.0, 1.0], [0.0, 1.0], [1.0, 0.0]]),
                         observation_ids=['f2_1539', 'f3_44', '295053'],
                         sample_ids=['f2', 'f3'])
        self.assertEqual(actual, expected)

    def test_underscores_in_sample_id(self):
        """ sample id with underscores is correctly processed
        """
        actual = parse_uc(uc_underscores_in_sample_id.split('\n'))
        expected = Table(np.array([[2.0, 1.0], [0.0, 1.0], [1.0, 0.0]]),
                         observation_ids=['_f_2__1539', 'f_3_44', '295053'],
                         sample_ids=['_f_2_', 'f_3'])
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    main()
