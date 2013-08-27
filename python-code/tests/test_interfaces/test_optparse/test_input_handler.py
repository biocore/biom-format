#!/usr/bin/env python

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

import os
from shutil import rmtree
from tempfile import mkdtemp
from biom.unit_test import TestCase, main
from biom.interfaces.optparse.input_handler import (load_biom_table,
        load_biom_table_with_file_contents, load_json_document, load_metadata)
from biom.parse import MetadataMap
from biom.table import SparseOTUTable

class InputHandlerTests(TestCase):
    def setUp(self):
        self.output_dir = mkdtemp()

        self.biom_fp = os.path.join(self.output_dir, 'test.biom')
        with open(self.biom_fp, 'w') as f:
            f.write(biom1)
            
        self.md_fp = os.path.join(self.output_dir, 'md.txt')
        with open(self.md_fp, 'w') as f:
            f.write(sample_md1)

    def tearDown(self):
        rmtree(self.output_dir)

    def test_load_biom_table(self):
        """Correctly parses and loads a BIOM table."""
        obs = load_biom_table(self.biom_fp)
        self.assertEqual(type(obs), SparseOTUTable)

    def test_load_biom_table_with_file_contents(self):
        """Correctly parses and loads a BIOM table, also returning the file."""
        obs = load_biom_table_with_file_contents(self.biom_fp)
        self.assertEqual(len(obs), 2)
        self.assertEqual(type(obs[0]), SparseOTUTable)
        self.assertEqual(type(obs[1]), file)
        obs[1].close()

    def test_load_json_document(self):
        """Correctly parses and loads a JSON document."""
        obs = load_json_document(self.biom_fp)
        self.assertEqual(type(obs), dict)
        self.assertEqual(obs['type'], 'OTU table')

    def test_load_metadata(self):
        """Correctly parses and loads a metadata map."""
        obs = load_metadata(self.md_fp)
        self.assertEqual(type(obs), MetadataMap)
        self.assertEqual(obs['t1'], {'foo': 't;b;c'})

        obs = load_metadata(None)
        self.assertTrue(obs is None)


biom1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "QIIME 1.6.0-dev","date": "2013-02-09T09:30:11.550590","matrix_type": "sparse","matrix_element_type": "int","shape": [14, 9],"data": [[0,0,20],[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],[2,4,1],[2,5,1],[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],[8,4,1],[8,6,2],[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],[11,6,1],[11,8,4],[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],[13,2,4]],"rows": [{"id": "295053", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "42684", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria"]}},{"id": "None11", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None10", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None7", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None6", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None5", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None4", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None3", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None2", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None1", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "879972", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None9", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None8", "metadata": {"taxonomy": ["k__Bacteria"]}}],"columns": [{"id": "f2", "metadata": null},{"id": "f1", "metadata": null},{"id": "f3", "metadata": null},{"id": "f4", "metadata": null},{"id": "p2", "metadata": null},{"id": "p1", "metadata": null},{"id": "t1", "metadata": null},{"id": "not16S.1", "metadata": null},{"id": "t2", "metadata": null}]}"""

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


if __name__ == '__main__':
    main()
