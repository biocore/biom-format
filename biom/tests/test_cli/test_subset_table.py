# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import os
import unittest

from biom.cli.table_subsetter import _subset_table
from biom.parse import parse_biom_table


class TestSubsetTable(unittest.TestCase):
    def setUp(self):
        """Set up data for use in unit tests."""
        self.biom_str1 = biom1

    def test_subset_samples(self):
        """Correctly subsets samples in a table."""
        obs = _subset_table(json_table_str=self.biom_str1, axis='sample',
                            ids=['f4', 'f2'], hdf5_biom=None)
        obs = parse_biom_table(list(obs[0]))
        self.assertEqual(len(obs.ids()), 2)
        self.assertEqual(len(obs.ids(axis='observation')), 14)
        self.assertTrue('f4' in obs.ids())
        self.assertTrue('f2' in obs.ids())

    def test_subset_observations(self):
        """Correctly subsets observations in a table."""
        obs = _subset_table(json_table_str=self.biom_str1, axis='observation',
                            ids=['None2', '879972'], hdf5_biom=None)
        obs = parse_biom_table(list(obs[0]))
        self.assertEqual(len(obs.ids()), 9)
        self.assertEqual(len(obs.ids(axis='observation')), 2)
        self.assertTrue('None2' in obs.ids(axis='observation'))
        self.assertTrue('879972' in obs.ids(axis='observation'))

    def test_invalid_input(self):
        """Correctly raises politically correct error upon invalid input."""
        with self.assertRaises(ValueError):
            _subset_table(hdf5_biom=None, json_table_str=self.biom_str1,
                          axis='foo', ids=['f2', 'f4'])

        with self.assertRaises(ValueError):
            _subset_table(hdf5_biom=None, json_table_str=None, axis='sample',
                          ids=['f2', 'f4'])

        with self.assertRaises(ValueError):
            _subset_table(json_table_str=self.biom_str1, hdf5_biom='foo',
                          axis='sample', ids=['f2', 'f4'])

    def test_subset_samples_hdf5(self):
        """Correctly subsets samples in a hdf5 table"""
        cwd = os.getcwd()
        if os.path.sep in __file__:
            os.chdir(os.path.dirname(__file__))
        obs = _subset_table(hdf5_biom=os.path.join('test_data', 'test.biom'),
                            axis='sample',
                            ids=['Sample1', 'Sample2', 'Sample3'],
                            json_table_str=None)
        os.chdir(cwd)
        obs = obs[0]
        self.assertEqual(len(obs.ids()), 3)
        self.assertEqual(len(obs.ids(axis='observation')), 5)
        self.assertTrue('Sample1' in obs.ids())
        self.assertTrue('Sample2' in obs.ids())
        self.assertTrue('Sample3' in obs.ids())

    def test_subset_observations_hdf5(self):
        """Correctly subsets samples in a hdf5 table"""
        cwd = os.getcwd()
        if os.path.sep in __file__:
            os.chdir(os.path.dirname(__file__))
        obs = _subset_table(hdf5_biom=os.path.join('test_data', 'test.biom'),
                            axis='observation',
                            ids=['GG_OTU_1', 'GG_OTU_3', 'GG_OTU_5'],
                            json_table_str=None)
        os.chdir(cwd)
        obs = obs[0]
        self.assertEqual(len(obs.ids()), 4)
        self.assertEqual(len(obs.ids(axis='observation')), 3)
        self.assertTrue('GG_OTU_1' in obs.ids(axis='observation'))
        self.assertTrue('GG_OTU_3' in obs.ids(axis='observation'))
        self.assertTrue('GG_OTU_5' in obs.ids(axis='observation'))


biom1 = ('{"id": "None","format": "Biological Observation Matrix 1.0.0",'
         '"format_url": "http://biom-format.org","type": "OTU table",'
         '"generated_by": "QIIME 1.6.0-dev","date": '
         '"2013-02-09T09:30:11.550590","matrix_type": "sparse",'
         '"matrix_element_type": "int","shape": [14, 9],"data": '
         '[[0,0,20],[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],[2,4,1]'
         ',[2,5,1],[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],[8,4,1],'
         '[8,6,2],[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],[11,6,1],'
         '[11,8,4],[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],[13,2,4]],'
         '"rows": [{"id": "295053", "metadata": {"taxonomy": ["k__Bacteria"]}}'
         ',{"id": "42684", "metadata": {"taxonomy": ["k__Bacteria", '
         '"p__Proteobacteria"]}},{"id": "None11", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "None10", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "None7", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "None6", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "None5", "metadata": {"taxonomy": '
         '["k__Bacteria"]}},{"id": "None4", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "None3", "metadata": {"taxonomy": '
         '["k__Bacteria"]}},{"id": "None2", "metadata": {"taxonomy": '
         '["k__Bacteria"]}},{"id": "None1", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "879972", "metadata": {"taxonomy": '
         '["k__Bacteria"]}},{"id": "None9", "metadata": {"taxonomy": '
         '["Unclassified"]}},{"id": "None8", "metadata": {"taxonomy": '
         '["k__Bacteria"]}}],"columns": [{"id": "f2", "metadata": null},'
         '{"id": "f1", "metadata": null},{"id": "f3", "metadata": null},'
         '{"id": "f4", "metadata": null},{"id": "p2", "metadata": null},{"id":'
         ' "p1", "metadata": null},{"id": "t1", "metadata": null},{"id": '
         '"not16S.1", "metadata": null},{"id": "t2", "metadata": null}]}')

if __name__ == "__main__":
    unittest.main()
