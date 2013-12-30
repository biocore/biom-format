#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from biom.commands.table_summarizer import TableSummarizer
from biom.parse import parse_biom_table
from biom.unit_test import TestCase, main

class TableSummarizerTests(TestCase):
    
    def setUp(self):
        """ initialize objects for use in tests """
        self.biom1_lines = biom1
        self.summary_default_lines = summary_default.split('\n')
        self.summary_suppress_md5_lines = summary_suppress_md5.split('\n')
        self.summary_qualitative_lines = summary_qualitative.split('\n')
    
    def test_default(self):
        """ TableSummarizer functions as expected with md5 suppression
        
        """
        t = TableSummarizer()
        actual = t(table=(parse_biom_table(self.biom1_lines),
                          self.biom1_lines.split('\n')),
                          qualitative=False,
                          suppress_md5=False)
        self.assertEqual(actual['biom_summary'],self.summary_default_lines)
    
    def test_suppress_md5(self):
        """ TableSummarizer functions as expected with md5 suppression
        
        """
        t = TableSummarizer()
        # suppress md5 by passing suppress_md5=True
        actual = t(table=(parse_biom_table(self.biom1_lines),
                                           self.biom1_lines),
                   qualitative=False,
                   suppress_md5=True)
        
        self.assertEqual(actual['biom_summary'],
                         self.summary_suppress_md5_lines)
        # suppress md5 by passing None as the second value in table
        actual = t(table=(parse_biom_table(self.biom1_lines),None),
                   qualitative=False,
                   suppress_md5=False)
        self.assertEqual(actual['biom_summary'],
                         self.summary_suppress_md5_lines)
        # suppress md5 by passing None as the second value in table 
        # and suppress_md5=True
        actual = t(table=(parse_biom_table(self.biom1_lines),None),
                   qualitative=False,
                   suppress_md5=True)
        self.assertEqual(actual['biom_summary'],
                         self.summary_suppress_md5_lines)

    def test_qualitative(self):
        """ TableSummarizer functions as expected with qualitative=True
        
        """
        t = TableSummarizer()
        actual = t(table=(parse_biom_table(self.biom1_lines),
                          self.biom1_lines.split('\n')),
                   qualitative=True,
                   suppress_md5=False)
        self.assertEqual(actual['biom_summary'],self.summary_qualitative_lines)

biom1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "QIIME 1.6.0-dev","date": "2013-02-09T09:30:11.550590","matrix_type": "sparse","matrix_element_type": "int","shape": [14, 9],"data": [[0,0,20],[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],[2,4,1],[2,5,1],[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],[8,4,1],[8,6,2],[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],[11,6,1],[11,8,4],[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],[13,2,4]],"rows": [{"id": "295053", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "42684", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria"]}},{"id": "None11", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None10", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None7", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None6", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None5", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None4", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None3", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None2", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None1", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "879972", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "None9", "metadata": {"taxonomy": ["Unclassified"]}},{"id": "None8", "metadata": {"taxonomy": ["k__Bacteria"]}}],"columns": [{"id": "f2", "metadata": null},{"id": "f1", "metadata": null},{"id": "f3", "metadata": null},{"id": "f4", "metadata": null},{"id": "p2", "metadata": null},{"id": "p1", "metadata": null},{"id": "t1", "metadata": null},{"id": "not16S.1", "metadata": null},{"id": "t2", "metadata": null}]}"""

summary_default = """Num samples: 9
Num observations: 14
Total count: 200
Table density (fraction of non-zero values): 0.238
Table md5 (unzipped): 27c6ffe253527068c82cc1835adeda8f

Counts/sample summary:
 Min: 22.0
 Max: 23.0
 Median: 22.000
 Mean: 22.222
 Std. dev.: 0.416
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
 p2: 22.0
 f1: 22.0
 f2: 22.0
 f3: 22.0
 f4: 22.0
 t2: 22.0
 not16S.1: 22.0
 t1: 23.0
 p1: 23.0"""

summary_suppress_md5 = """Num samples: 9
Num observations: 14
Total count: 200
Table density (fraction of non-zero values): 0.238

Counts/sample summary:
 Min: 22.0
 Max: 23.0
 Median: 22.000
 Mean: 22.222
 Std. dev.: 0.416
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
 p2: 22.0
 f1: 22.0
 f2: 22.0
 f3: 22.0
 f4: 22.0
 t2: 22.0
 not16S.1: 22.0
 t1: 23.0
 p1: 23.0"""

summary_qualitative = """Num samples: 9
Num observations: 14
Table md5 (unzipped): 27c6ffe253527068c82cc1835adeda8f

Observations/sample summary:
 Min: 1
 Max: 9
 Median: 3.000
 Mean: 3.333
 Std. dev.: 2.211
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Observations/sample detail:
 f4: 1
 f1: 2
 f3: 2
 not16S.1: 2
 f2: 3
 t2: 3
 t1: 4
 p1: 4
 p2: 9"""

if __name__ == "__main__":
    main()
