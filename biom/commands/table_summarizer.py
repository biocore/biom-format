#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from pyqi.core.command import (Command, CommandIn, CommandOut, 
        ParameterCollection)

from numpy import std
from operator import itemgetter
from biom.util import (compute_counts_per_sample_stats, 
                       safe_md5,
                       biom_open)
from biom.parse import parse_biom_table
from biom.table import Table

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Daniel McDonald"]
__license__ = "BSD"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class TableSummarizer(Command):
    """
     Example usage:
      from biom.commands.table_summarizer import TableSummarizer
      from biom.parse import parse_biom_table
      c = TableSummarizer()
      table_f = open("table.biom")
      t = parse_biom_table(table_f)
      table_f.seek(0)
      result = c(table=(t,None))
      result = c(table=(t,None),qualitative=True)
      result = c(table=(t,table_f),qualitative=True)
      table_f.close()
    """
    BriefDescription = "Summarize sample or observation data in a BIOM table"
    LongDescription = "Provides details on the observation counts per sample, "\
                      "including summary statistics, as well as metadata "\
                      "categories associated with samples and observations."
    
    CommandIns = ParameterCollection([
        CommandIn(Name='table', 
                  DataType=tuple,
                  Description='the input BIOM table', 
                  Required=True),
        CommandIn(Name='qualitative', 
                  DataType=bool,
                  Description=('Present counts as number of unique '
                               'observation ids per sample, rather than '
                               'counts of observations per sample.'), 
                  Required=False,
                  Default=False),
        CommandIn(Name='suppress_md5', 
                  DataType=bool,
                  Description=('Do not compute md5sum of table. '
                               'Useful if you\'re concerned about runtime.'), 
                  Required=False,
                  Default=False)
    ])
    
    CommandOuts = ParameterCollection([
        CommandOut(Name='biom_summary',
                   DataType=list,
                   Description='The table summary')
    ])
                  
    def run(self, **kwargs):
        result = {}
        qualitative = kwargs['qualitative']
        table, table_lines = kwargs['table']
       
        min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
            compute_counts_per_sample_stats(table, qualitative)
        num_observations = len(table.ObservationIds)
        
        suppress_md5 = (table_lines is None) or kwargs['suppress_md5']
    
        counts_per_sample_values = counts_per_sample.values()
    
        if table.SampleMetadata is None:
            sample_md_keys = ["None provided"]
        else:
            sample_md_keys = table.SampleMetadata[0].keys()
        
        if table.ObservationMetadata is None:
            observation_md_keys = ["None provided"]
        else:
            observation_md_keys = table.ObservationMetadata[0].keys()
    
        lines = []
    
        num_samples = len(table.SampleIds)
        lines.append('Num samples: %d' % num_samples)
        lines.append('Num observations: %d' % num_observations)
        
        if not qualitative:
            total_count = sum(counts_per_sample_values)
            lines.append('Total count: %d' % total_count)
            lines.append('Table density (fraction of non-zero values): %1.3f' % \
                  table.getTableDensity())
        
        if not suppress_md5:
            lines.append('Table md5 (unzipped): %s' % safe_md5(table_lines))
        lines.append('')

        if qualitative:
            lines.append('Observations/sample summary:')
        else:
            lines.append('Counts/sample summary:')
        
        lines.append(' Min: %r' % min_counts)
        lines.append(' Max: %r' % max_counts)
        lines.append(' Median: %1.3f' % median_counts)
        lines.append(' Mean: %1.3f' % mean_counts)
        lines.append(' Std. dev.: %1.3f' % std(counts_per_sample_values))
        lines.append(' Sample Metadata Categories: %s' % '; '.join(sample_md_keys))
        lines.append(' Observation Metadata Categories: %s' % '; '.join(observation_md_keys))
        lines.append('')
        
        if qualitative:
            lines.append('Observations/sample detail:')
        else:
            lines.append('Counts/sample detail:')
        
        for k,v in sorted(counts_per_sample.items(), key=itemgetter(1)):
            lines.append(' %s: %r' % (k,v))
        
        result['biom_summary'] = lines
        return result
    
CommandConstructor = TableSummarizer
