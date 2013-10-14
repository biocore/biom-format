#!/usr/bin/env python

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from pyqi.core.exception import CommandError
from biom.commands.table_subsetter import TableSubsetter
from biom.parse import parse_biom_table
from biom.unit_test import TestCase, main

class TableSubsetterTests(TestCase):
    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = TableSubsetter()
        self.biom_str1 = biom1

    def test_subset_samples(self):
        """Correctly subsets samples in a table."""
        obs = self.cmd(table_str=self.biom_str1, axis='samples',
                       ids=['f4', 'f2'])
        self.assertEqual(obs.keys(), ['subset_generator'])
        obs = parse_biom_table(list(obs['subset_generator']))
        self.assertEqual(len(obs.SampleIds), 2)
        self.assertEqual(len(obs.ObservationIds), 14)
        self.assertTrue('f4' in obs.SampleIds)
        self.assertTrue('f2' in obs.SampleIds)

    def test_subset_observations(self):
        """Correctly subsets observations in a table."""
        obs = self.cmd(table_str=self.biom_str1, axis='observations',
                       ids=['None2', '879972'])
        self.assertEqual(obs.keys(), ['subset_generator'])
        obs = parse_biom_table(list(obs['subset_generator']))
        self.assertEqual(len(obs.SampleIds), 9)
        self.assertEqual(len(obs.ObservationIds), 2)
        self.assertTrue('None2' in obs.ObservationIds)
        self.assertTrue('879972' in obs.ObservationIds)

    def test_invalid_input(self):
        """Correctly raises politically correct error upon invalid input."""
        with self.assertRaises(CommandError):
            _ = self.cmd(table_str=self.biom_str1, axis='foo',
                         ids=['f2', 'f4'])


biom1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "QIIME 1.6.0-dev","date": "2013-02-09T09:30:11.550590","matrix_type": "sparse","matrix_element_type": "int","shape": [14, 9],"data": [[0,0,20],[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],[2,4,1],[2,5,1],[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],[8,4,1],[8,6,2],[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],[11,6,1],[11,8,4],[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],[13,2,4]],"rows": [{"id": "295053", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "42684", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria"]}},{"id": "None11", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None10", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None7", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None6", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None5", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None4", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None3", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None2", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None1", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "879972", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None9", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None8", "metadata": {"taxonomy": ["k__Bacteria"]}}],"columns": [{"id": "f2", "metadata": null},{"id": "f1", "metadata": null},{"id": "f3", "metadata": null},{"id": "f4", "metadata": null},{"id": "p2", "metadata": null},{"id": "p1", "metadata": null},{"id": "t1", "metadata": null},{"id": "not16S.1", "metadata": null},{"id": "t2", "metadata": null}]}"""


if __name__ == "__main__":
    main()
