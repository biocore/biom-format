#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from pyqi.core.exception import CommandError
from biom.commands.metadata_adder import MetadataAdder
from biom.parse import parse_biom_table
from biom.unit_test import TestCase, main

class MetadataAdderTests(TestCase):
    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = MetadataAdder()
        self.biom_lines1 = biom1
        self.biom_table1 = parse_biom_table(self.biom_lines1)
        self.sample_md_lines1 = sample_md1.split('\n')
        self.obs_md_lines1 = obs_md1.split('\n')

    def test_add_sample_metadata_no_casting(self):
        """Correctly adds sample metadata without casting it."""
        # Add a subset of sample metadata to a table that doesn't have any
        # sample metadata to begin with. Don't perform any casting.
        obs = self.cmd(table=self.biom_table1,
                       sample_metadata=self.sample_md_lines1)
        self.assertEqual(obs.keys(), ['table'])

        obs = obs['table']
        self.assertEqual(obs.SampleMetadata[obs.getSampleIndex('f4')],
                         {'bar': '0.23', 'foo': '9', 'baz': 'abc;123'})
        self.assertEqual(obs.SampleMetadata[obs.getSampleIndex('not16S.1')],
                         {'bar': '-4.2', 'foo': '0', 'baz': '123;abc'})
        self.assertEqual(obs.SampleMetadata[obs.getSampleIndex('f2')], {})

    def test_add_sample_metadata_with_casting(self):
        """Correctly adds sample metadata with casting."""
        obs = self.cmd(table=self.biom_table1,
                       sample_metadata=self.sample_md_lines1,
                       sc_separated=['baz'], int_fields=['foo'],
                       float_fields=['bar'])
        self.assertEqual(obs.keys(), ['table'])

        obs = obs['table']
        self.assertEqual(obs.SampleMetadata[obs.getSampleIndex('f4')],
                         {'bar': 0.23, 'foo': 9, 'baz': ['abc', '123']})
        self.assertEqual(obs.SampleMetadata[obs.getSampleIndex('not16S.1')],
                         {'bar': -4.2, 'foo': 0, 'baz': ['123', 'abc']})
        self.assertEqual(obs.SampleMetadata[obs.getSampleIndex('f2')], {})

    def test_add_observation_metadata_no_casting(self):
        """Correctly adds observation metadata without casting it."""
        # Add observation metadata to a table that already has observation
        # metadata. Some observations won't be modified, and metadata for
        # observations that aren't in the table are included. Don't perform any
        # casting.
        obs = self.cmd(table=self.biom_table1,
                       observation_metadata=self.obs_md_lines1)
        self.assertEqual(obs.keys(), ['table'])

        obs = obs['table']
        self.assertEqual(
                obs.ObservationMetadata[obs.getObservationIndex('None7')],
                         {'foo': '6', 'taxonomy': 'abc;123|def;456'})
        self.assertEqual(
                obs.ObservationMetadata[obs.getObservationIndex('879972')],
                         {'foo': '3', 'taxonomy': '123;abc|456;def'})
        self.assertEqual(
                obs.ObservationMetadata[obs.getObservationIndex('None8')],
                {'taxonomy': ['k__Bacteria']})

    def test_add_observation_metadata_with_casting(self):
        """Correctly adds observation metadata with casting."""
        obs = self.cmd(table=self.biom_table1,
                       observation_metadata=self.obs_md_lines1,
                       sc_pipe_separated=['taxonomy'], int_fields=['foo'])
        self.assertEqual(obs.keys(), ['table'])

        obs = obs['table']
        self.assertEqual(
                obs.ObservationMetadata[obs.getObservationIndex('None7')],
                         {'foo': 6, 'taxonomy': [['abc','123'],['def','456']]})
        self.assertEqual(
                obs.ObservationMetadata[obs.getObservationIndex('879972')],
                         {'foo': 3, 'taxonomy': [['123','abc'],['456','def']]})
        self.assertEqual(
                obs.ObservationMetadata[obs.getObservationIndex('None8')],
                {'taxonomy': ['k__Bacteria']})

    def test_no_metadata(self):
        """Correctly raises error when not provided any metadata."""
        with self.assertRaises(CommandError):
            _ = self.cmd(table=self.biom_table1)


biom1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "QIIME 1.6.0-dev","date": "2013-02-09T09:30:11.550590","matrix_type": "sparse","matrix_element_type": "int","shape": [14, 9],"data": [[0,0,20],[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],[2,4,1],[2,5,1],[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],[8,4,1],[8,6,2],[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],[11,6,1],[11,8,4],[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],[13,2,4]],"rows": [{"id": "295053", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "42684", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria"]}},{"id": "None11", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None10", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None7", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None6", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None5", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None4", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None3", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None2", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None1", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "879972", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None9", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None8", "metadata": {"taxonomy": ["k__Bacteria"]}}],"columns": [{"id": "f2", "metadata": null},{"id": "f1", "metadata": null},{"id": "f3", "metadata": null},{"id": "f4", "metadata": null},{"id": "p2", "metadata": null},{"id": "p1", "metadata": null},{"id": "t1", "metadata": null},{"id": "not16S.1", "metadata": null},{"id": "t2", "metadata": null}]}"""

sample_md1 = """#SampleID\tfoo\tbar\tbaz
f4\t9\t0.23\tabc;123
not16S.1\t0\t-4.2\t123;abc
"""

obs_md1 = """#OTUID\tfoo\ttaxonomy
None7\t6\tabc;123|def;456
best-observation\t8\tghi;789|jkl;101112
879972\t3\t123;abc|456;def
"""


if __name__ == "__main__":
    main()
