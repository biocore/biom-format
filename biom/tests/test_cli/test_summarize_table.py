#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from biom.cli.table_summarizer import _summarize_table
from biom.parse import load_table

import tempfile
import os
from unittest import TestCase, main


class TestSummarizeTable(TestCase):

    def setUp(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as fh:
            fh.write(biom1)
            fh.flush()
            self.biom1 = load_table(fh.name)
            self.temporary_fh_name = fh.name

    def tearDown(self):
        os.unlink(self.temporary_fh_name)

    def test_default(self):
        """ TableSummarizer functions as expected

        """
        result = _summarize_table(self.biom1)
        # test same alphanumeric content, order of samples is runtime
        # dependent
        self.assertEqual(sorted(result), sorted(summary_default))

    def test_qualitative(self):
        """ TableSummarizer functions as expected with qualitative=True

        """
        result = _summarize_table(self.biom1, qualitative=True)
        # test same alphanumeric content, order of samples is runtime
        # dependent
        self.assertEqual(sorted(result), sorted(summary_qualitative))


biom1 = ('{"id": "None","format": "Biological Observation Matrix 1.0.0",'
         '"format_url": "http://biom-format.org","type": "OTU table",'
         '"generated_by": "QIIME 1.6.0-dev","date": '
         '"2013-02-09T09:30:11.550590","matrix_type": "sparse",'
         '"matrix_element_type": "int","shape": [14, 9],"data": [[0,0,20],'
         '[0,1,18],[0,2,18],[0,3,22],[0,4,4],[1,4,1],[2,0,1],[2,4,1],[2,5,1],'
         '[3,6,1],[4,4,1],[5,7,20],[6,4,1],[7,4,1],[7,5,1],[8,4,1],[8,6,2],'
         '[8,8,3],[9,7,2],[10,5,1],[11,4,9],[11,5,20],[11,6,1],[11,8,4],'
         '[12,4,3],[12,6,19],[12,8,15],[13,0,1],[13,1,4],[13,2,4]],"rows": '
         '[{"id": "295053", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": '
         '"42684", "metadata": {"taxonomy": ["k__Bacteria", '
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

summary_default = """Num samples: 9
Num observations: 14
Total count: 200
Table density (fraction of non-zero values): 0.238

Counts/sample summary:
 Min: 22.000
 Max: 23.000
 Median: 22.000
 Mean: 22.222
 Std. dev.: 0.416
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
p2: 22.000
f1: 22.000
f2: 22.000
f3: 22.000
f4: 22.000
t2: 22.000
not16S.1: 22.000
t1: 23.000
p1: 23.000"""

summary_qualitative = """Num samples: 9
Num observations: 14

Observations/sample summary:
 Min: 1.000
 Max: 9.000
 Median: 3.000
 Mean: 3.333
 Std. dev.: 2.211
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Observations/sample detail:
f4: 1.000
f1: 2.000
f3: 2.000
not16S.1: 2.000
f2: 3.000
t2: 3.000
t1: 4.000
p1: 4.000
p2: 9.000"""

if __name__ == "__main__":
    main()
