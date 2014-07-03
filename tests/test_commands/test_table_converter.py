#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from pyqi.core.exception import CommandError
from biom.commands.table_converter import TableConverter
from biom.parse import MetadataMap, parse_biom_table
from biom.table import Table
from unittest import TestCase, main
from StringIO import StringIO


class TableConverterTests(TestCase):

    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = TableConverter()

        self.biom_lines1 = biom1
        self.biom_table1 = parse_biom_table(self.biom_lines1)

        self.classic_lines1 = classic1.split('\n')

        self.sample_md1 = MetadataMap.from_file(sample_md1.split('\n'))

    def test_classic_to_biom(self):
        """Correctly converts classic to biom."""
        obs = self.cmd(table=parse_biom_table(self.classic_lines1),
                       to_json=True, table_type='OTU table')
        self.assertEqual(obs.keys(), ['table'])

        obs = parse_biom_table(obs['table'][0].to_json('testing'))
        self.assertEqual(type(obs), Table)
        self.assertEqual(len(obs.ids()), 9)
        self.assertEqual(len(obs.ids(axis='observation')), 14)
        self.assertEqual(obs.metadata(), None)
        self.assertNotEqual(obs.metadata(axis='observation'), None)

    def test_classic_to_biom_with_metadata(self):
        """Correctly converts classic to biom with metadata."""
        # No processing of metadata.
        obs = self.cmd(table=parse_biom_table(self.classic_lines1),
                       sample_metadata=self.sample_md1, to_json=True,
                       table_type='OTU table', process_obs_metadata='naive')
        self.assertEqual(obs.keys(), ['table'])

        obs = parse_biom_table(obs['table'][0].to_json('testing'))
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
        obs = self.cmd(table=parse_biom_table(self.classic_lines1),
                       sample_metadata=self.sample_md1, table_type='OTU table',
                       process_obs_metadata='sc_separated', to_json=True)
        self.assertEqual(obs.keys(), ['table'])

        obs = parse_biom_table(obs['table'][0].to_json('testing'))
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

    def test_biom_to_classic(self):
        """Correctly converts biom to classic."""
        obs = self.cmd(table=parse_biom_table(self.biom_lines1),
                       to_tsv=True, header_key='taxonomy')
        self.assertEqual(obs.keys(), ['table'])
        self.assertEqual(obs['table'][0], classic1)

        obs = self.cmd(table=parse_biom_table(self.biom_lines1), to_tsv=True,
                       header_key='taxonomy', output_metadata_id='foo')
        self.assertEqual(obs.keys(), ['table'])
        obs_md_col = obs['table'][0].split('\n')[1].split('\t')[-1]
        self.assertEqual(obs_md_col, 'foo')

    def test_invalid_input(self):
        """Correctly handles invalid input by raising a CommandError."""
        # Too many ops.
        with self.assertRaises(CommandError):
            self.cmd(table_file=self.biom_lines1,
                     sparse_biom_to_dense_biom=True,
                     biom_to_classic_table=True)

        # biom -> classic, but supply classic
        with self.assertRaises(CommandError):
            self.cmd(table_file=self.classic_lines1,
                     biom_to_classic_table=True)

        # sparse biom -> dense biom, but supply classic
        with self.assertRaises(CommandError):
            self.cmd(table_file=self.classic_lines1,
                     sparse_biom_to_dense_biom=True)

        # dense biom -> sparse biom, but supply classic
        with self.assertRaises(CommandError):
            self.cmd(table_file=self.classic_lines1,
                     dense_biom_to_sparse_biom=True)

        # Unknown observation processor.
        with self.assertRaises(CommandError):
            self.cmd(table_file=self.classic_lines1,
                     process_obs_metadata='foo')

        # classic -> biom, but supply biom
        with self.assertRaises(CommandError):
            self.cmd(table_file=StringIO(self.biom_lines1),
                     process_obs_metadata='sc_separated')


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
