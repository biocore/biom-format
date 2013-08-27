#!/usr/bin/env python

from __future__ import division
from pyqi.core.command import Command, Parameter, ParameterCollection

from numpy import std

from biom.util import (compute_counts_per_sample_stats, 
                       safe_md5,
                       biom_open)
from biom.parse import parse_biom_table
from biom.table import Table

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class TableSummarizer(Command):
    """
     Example usage:
      from biom.commands.table_summarizer import TableSummarizer
      from biom.parse import parse_biom_table
      c = TableSummarizer()
      table_f = open("table.biom","U")
      t = parse_biom_table(table_f)
      table_f.seek(0)
      result = c(table=(t,None))
      result = c(table=(t,None),qualitative=True)
      result = c(table=(t,table_f),qualitative=True)
      table_f.close()
    """
    BriefDescription = "Summarize sample or observation data in a BIOM table"
    LongDescription = "Provides details on the observation counts per sample, including summary statistics, as well as metadata categories associated with samples and observations."
    Parameters = ParameterCollection([
        Parameter(Name='table', 
                  DataType=tuple,
                  Description='the input BIOM table', 
                  Required=True),
        Parameter(Name='qualitative', 
                  DataType=bool,
                  Description=('Present counts as number of unique '
                               'observation ids per sample, rather than '
                               'counts of observations per sample.'), 
                  Required=False,
                  Default=False),
        Parameter(Name='suppress_md5', 
                  DataType=bool,
                  Description=('Do not compute md5sum of table. '
                               'Useful if you\'re concerned about runtime.'), 
                  Required=False,
                  Default=False)
    ])

    def run(self, **kwargs):
        """
         table: two-element tuple containing the biom table to summarize and
                the file(-like) object containing the original table data. The
                second element of the tuple (the file(-like) object) may be
                None. If this is the case, the MD5 sum will *not* be computed
         qualitative: counts are presented as number of unique observation
                      ids per sample, rather than total observation count per
                      sample
         suppress_md5: if ``True``, the MD5 sum of the table file contents will
                       not be computed. This parameter is ignored if
                       ``table[1] is None``
        """
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
    
        num_samples = len(counts_per_sample)
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
        
        sorted_counts_per_sample = [(v,k) for k,v in counts_per_sample.items()]
        sorted_counts_per_sample.sort()
        for v,k in sorted_counts_per_sample:
            lines.append(' %s: %r' % (k,v))
        
        result['biom-summary'] = lines
        return result
    

CommandConstructor = TableSummarizer
