#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import os
from os.path import abspath, dirname, join
import tempfile

import numpy as np

from biom.cli.table_converter import _convert
from biom.parse import MetadataMap
from biom.table import Table
from biom import load_table
from biom.parse import biom_open, parse_biom_table
from unittest import TestCase, main


class TableConverterTests(TestCase):

    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = _convert
        self.output_filepath = tempfile.NamedTemporaryFile().name

        with tempfile.NamedTemporaryFile('w', delete=False) as fh:
            fh.write(biom1)
            fh.flush()
            self.biom_table1 = load_table(fh.name)
            self.temporary_fh_table_name = fh.name

        self.biom_lines1 = biom1.split('\n')
        with tempfile.NamedTemporaryFile('w', delete=False) as fh:
            fh.write(classic1)
            fh.flush()
            self.classic_biom1 = load_table(fh.name)
            self.temporary_fh_classic_name = fh.name

        self.sample_md1 = MetadataMap.from_file(sample_md1.split('\n'))

        test_data_dir = join(dirname(abspath(__file__)), 'test_data')
        self.json_collapsed_obs = join(test_data_dir,
                                       'json_obs_collapsed.biom')
        self.json_collapsed_samples = join(test_data_dir,
                                           'json_sample_collapsed.biom')

    def tearDown(self):
        os.unlink(self.temporary_fh_classic_name)
        os.unlink(self.temporary_fh_table_name)

    def test_classic_to_biom(self):
        """Correctly converts classic to biom."""
        self.cmd(table=self.classic_biom1,
                 output_filepath=self.output_filepath,
                 to_json=True, table_type='OTU table')

        obs = load_table(self.output_filepath)
        self.assertEqual(type(obs), Table)
        self.assertEqual(len(obs.ids()), 9)
        self.assertEqual(len(obs.ids(axis='observation')), 14)
        self.assertEqual(obs.metadata(), None)
        self.assertNotEqual(obs.metadata(axis='observation'), None)

    def test_classic_to_biom_with_metadata(self):
        """Correctly converts classic to biom with metadata."""
        # No processing of metadata.
        obs = self.cmd(table=self.classic_biom1,
                       output_filepath=self.output_filepath,
                       sample_metadata=self.sample_md1, to_json=True,
                       table_type='OTU table', process_obs_metadata='naive')

        obs = load_table(self.output_filepath)
        self.assertEqual(type(obs), Table)
        self.assertEqual(len(obs.ids()), 9)
        self.assertEqual(len(obs.ids(axis='observation')), 14)
        self.assertNotEqual(obs.metadata(), None)
        self.assertNotEqual(obs.metadata(axis='observation'), None)
        self.assertEqual(obs.metadata()[obs.index('p2', 'sample')],
                         {'foo': 'c;b;a'})
        self.assertEqual(obs.metadata()[obs.index('not16S.1', 'sample')],
                         {'foo': 'b;c;d'})
        self.assertEqual(obs.metadata(axis='observation')[
            obs.index('None11', 'observation')],
            {'taxonomy': 'Unclassified'})

        # With processing of metadata (currently only supports observation md).
        obs = self.cmd(table=self.classic_biom1,
                       output_filepath=self.output_filepath,
                       sample_metadata=self.sample_md1, table_type='OTU table',
                       process_obs_metadata='sc_separated', to_json=True)

        obs = load_table(self.output_filepath)
        self.assertEqual(type(obs), Table)
        self.assertEqual(len(obs.ids()), 9)
        self.assertEqual(len(obs.ids(axis='observation')), 14)
        self.assertNotEqual(obs.metadata(), None)
        self.assertNotEqual(obs.metadata(axis='observation'), None)
        self.assertEqual(obs.metadata()[obs.index('p2', 'sample')],
                         {'foo': 'c;b;a'})
        self.assertEqual(obs.metadata()[obs.index('not16S.1', 'sample')],
                         {'foo': 'b;c;d'})
        self.assertEqual(obs.metadata(axis='observation')[
            obs.index('None11', 'observation')],
            {'taxonomy': ['Unclassified']})

    def test_biom_to_classic1(self):
        """Correctly converts biom to classic."""
        self.cmd(table=self.biom_table1,
                 output_filepath=self.output_filepath,
                 to_tsv=True, header_key='taxonomy')

        self.assertEqual(load_table(self.output_filepath), self.classic_biom1)

    def test_biom_to_classic2(self):
        """Correctly converts biom to classic with metadata renaming."""
        self.cmd(table=self.biom_table1,
                 output_filepath=self.output_filepath, to_tsv=True,
                 header_key='taxonomy', output_metadata_id='foo')
        obs = load_table(self.output_filepath)
        self.assertTrue('foo' in obs.metadata(axis='observation')[0])

    def test_json_to_hdf5_collapsed_samples(self):
        """Correctly converts json to HDF5 changing the sample metadata"""
        with biom_open(self.json_collapsed_samples) as f:
            obs = self.cmd(table=parse_biom_table(f),
                           output_filepath=self.output_filepath, to_hdf5=True,
                           collapsed_samples=True)
        obs = load_table(self.output_filepath)
        exp = Table(np.array([[0., 1.], [6., 6.], [6., 1.],
                              [1., 4.], [0., 2.]]),
                    observation_ids=[
                        'GG_OTU_1',
                        'GG_OTU_2',
                        'GG_OTU_3',
                        'GG_OTU_4',
                        'GG_OTU_5',
                    ],
                    sample_ids=['skin', 'gut'],
                    observation_metadata=[
                        {'taxonomy': [
                            'k__Bacteria',
                            'p__Proteobacteria',
                            'c__Gammaproteobacteria',
                            'o__Enterobacteriales',
                            'f__Enterobacteriaceae',
                            'g__Escherichia',
                            's__',
                        ]},
                        {'taxonomy': [
                            'k__Bacteria',
                            'p__Cyanobacteria',
                            'c__Nostocophycideae',
                            'o__Nostocales',
                            'f__Nostocaceae',
                            'g__Dolichospermum',
                            's__',
                        ]},
                        {'taxonomy': [
                            'k__Archaea',
                            'p__Euryarchaeota',
                            'c__Methanomicrobia',
                            'o__Methanosarcinales',
                            'f__Methanosarcinaceae',
                            'g__Methanosarcina',
                            's__',
                        ]},
                        {'taxonomy': [
                            'k__Bacteria',
                            'p__Firmicutes',
                            'c__Clostridia',
                            'o__Halanaerobiales',
                            'f__Halanaerobiaceae',
                            'g__Halanaerobium',
                            's__Halanaerobiumsaccharolyticum',
                        ]},
                        {'taxonomy': [
                            'k__Bacteria',
                            'p__Proteobacteria',
                            'c__Gammaproteobacteria',
                            'o__Enterobacteriales',
                            'f__Enterobacteriaceae',
                            'g__Escherichia',
                            's__',
                        ]}],
                    sample_metadata=[
                        {'collapsed_ids': [
                            'Sample4',
                            'Sample5',
                            'Sample6',
                        ]},
                        {'collapsed_ids': [
                            'Sample1',
                            'Sample2',
                            'Sample3',
                        ]}
                    ],
                    type='OTU table')
        self.assertEqual(obs, exp)

    def test_json_to_hdf5_collapsed_metadata(self):
        """Correctly converts json to HDF5 changing the observation metadata"""
        with biom_open(self.json_collapsed_obs) as f:
            t = parse_biom_table(f)
            obs = self.cmd(table=t,
                           output_filepath=self.output_filepath, to_hdf5=True,
                           collapsed_observations=True)
        obs = load_table(self.output_filepath)
        exp = Table(np.array([[2., 1., 1., 0., 0., 1.],
                              [0., 0., 1., 4., 0., 2.],
                              [5., 1., 0., 2., 3., 1.],
                              [0., 1., 2., 0., 0., 0.]]),
                    observation_ids=[
                        'p__Firmicutes',
                        'p__Euryarchaeota',
                        'p__Cyanobacteria',
                        'p__Proteobacteria',
                    ],
                    sample_ids=[
                        'Sample1',
                        'Sample2',
                        'Sample3',
                        'Sample4',
                        'Sample5',
                        'Sample6',
                    ],
                    observation_metadata=[
                        {'collapsed_ids': ['GG_OTU_4']},
                        {'collapsed_ids': ['GG_OTU_3']},
                        {'collapsed_ids': ['GG_OTU_2']},
                        {'collapsed_ids': ['GG_OTU_1', 'GG_OTU_5']}],
                    sample_metadata=[
                        {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                         'BarcodeSequence': 'CGCTTATCGAGA',
                         'Description': 'human gut',
                         'BODY_SITE': 'gut'},
                        {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                         'BarcodeSequence': 'CATACCAGTAGC',
                         'Description': 'human gut',
                         'BODY_SITE': 'gut'},
                        {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                         'BarcodeSequence': 'CTCTCTACCTGT',
                         'Description': 'human gut',
                         'BODY_SITE': 'gut'},
                        {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                         'BarcodeSequence': 'CTCTCGGCCTGT',
                         'Description': 'human skin',
                         'BODY_SITE': 'skin'},
                        {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                         'BarcodeSequence': 'CTCTCTACCAAT',
                         'Description': 'human skin',
                         'BODY_SITE': 'skin'},
                        {'LinkerPrimerSequence': 'CATGCTGCCTCCCGTAGGAGT',
                         'BarcodeSequence': 'CTAACTACCAAT',
                         'Description': 'human skin',
                         'BODY_SITE': 'skin'}],
                    type='OTU table')

        self.assertEqual(obs, exp)


biom1 = """
{"id": "None",
 "format": "Biological Observation Matrix 1.0.0",
 "format_url": "http://biom-format.org",
 "type": "OTU table",
 "generated_by": "QIIME 1.6.0-dev",
 "date": "2013-02-09T09:30:11.550590",
 "matrix_type": "sparse",
 "matrix_element_type": "float",
 "shape": [14, 9],
 "data": [[0,0,20],[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],
         [2,4,1],[2,5,1],[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],
         [8,4,1],[8,6,2],[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],
         [11,6,1],[11,8,4],[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],
         [13,2,4]],
 "rows": [{"id": "295053",
           "metadata": {"taxonomy": ["k__Bacteria"]}},
          {"id": "42684", "metadata": {"taxonomy": ["k__Bacteria",
                                                    "p__Proteobacteria"]}},
          {"id": "None11", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "None10", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "None7", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "None6", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "None5", "metadata": {"taxonomy": ["k__Bacteria"]}},
          {"id": "None4", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "None3", "metadata": {"taxonomy": ["k__Bacteria"]}},
          {"id": "None2", "metadata": {"taxonomy": ["k__Bacteria"]}},
          {"id": "None1", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "879972", "metadata": {"taxonomy": ["k__Bacteria"]}},
          {"id": "None9", "metadata": {"taxonomy": ["Unclassified"]}},
          {"id": "None8", "metadata": {"taxonomy": ["k__Bacteria"]}}],
 "columns": [{"id": "f2", "metadata": null},
             {"id": "f1", "metadata": null},{"id": "f3", "metadata": null},
             {"id": "f4", "metadata": null},{"id": "p2", "metadata": null},
             {"id": "p1", "metadata": null},{"id": "t1", "metadata": null},
             {"id": "not16S.1", "metadata": null},
             {"id": "t2", "metadata": null}]
 }"""

classic1 = """# Constructed from biom file
#OTU ID\tf2\tf1\tf3\tf4\tp2\tp1\tt1\tnot16S.1\tt2\ttaxonomy
295053\t20.0\t18.0\t18.0\t22.0\t4.0\t0.0\t0.0\t0.0\t0.0\tk__Bacteria
42684\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t0.0\tk__Bacteria; """ + \
    """p__Proteobacteria
None11\t1.0\t0.0\t0.0\t0.0\t1.0\t1.0\t0.0\t0.0\t0.0\tUnclassified
None10\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\tUnclassified
None7\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t0.0\tUnclassified
None6\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t20.0\t0.0\tUnclassified
None5\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t0.0\tk__Bacteria
None4\t0.0\t0.0\t0.0\t0.0\t1.0\t1.0\t0.0\t0.0\t0.0\tUnclassified
None3\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t2.0\t0.0\t3.0\tk__Bacteria
None2\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t2.0\t0.0\tk__Bacteria
None1\t0.0\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\tUnclassified
879972\t0.0\t0.0\t0.0\t0.0\t9.0\t20.0\t1.0\t0.0\t4.0\tk__Bacteria
None9\t0.0\t0.0\t0.0\t0.0\t3.0\t0.0\t19.0\t0.0\t15.0\tUnclassified
None8\t1.0\t4.0\t4.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\tk__Bacteria"""

sample_md1 = """#SampleID\tfoo
f4\ta;b;c
not16S.1\tb;c;d
f2\ta;c;d
f1\ta;b;c
p2\tc;b;a
f3\ta;b;c
t1\tt;b;c
p1\tp;b;c
t2\tt;2;z
"""


if __name__ == "__main__":
    main()
