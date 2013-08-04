#!/usr/bin/env python

from __future__ import division
from pyqi.core.command import Command, Parameter
from numpy import std

from biom.util import (compute_counts_per_sample_stats, 
                       safe_md5,
                       biom_open)
from biom.parse import parse_biom_table
from biom.table import Table

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Greg Caporaso", ""]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class PrintBiomTableSummary(Command):
    BriefDescription = "Summarize sample or observation data in a BIOM table."
    LongDescription = "Provides details on the observation counts per sample, including summary statistics, as well as metadata categories associated with samples and observations."

    def run(self, **kwargs):
        """
         table: the biom table to summarize
         input_fp: path to .biom file that table is derived from, provide if 
                   the md5 sum should be included in the output
         qualitative: counts are presented as number of unique observation
                   ids per sample, rather than total observation count per
                   sample
        """
        print kwargs
        result = {}
        qualitative = kwargs.get('qualitative')
        table = kwargs.get('table')
        
        min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
         compute_counts_per_sample_stats(table, qualitative)
        num_observations = len(table.ObservationIds)
        
        input_fp = kwargs.get('input_fp',None)
        suppress_md5 = input_fp is None
    
        counts_per_sample_values = counts_per_sample.values()
    
        try:
            sample_md_keys = table.SampleMetadata[0].keys()
        except TypeError:
            sample_md_keys = ["None provided"]
        try:
            observation_md_keys = table.ObservationMetadata[0].keys()
        except TypeError:
            observation_md_keys = ["None provided"]
    
        lines = []
    
        num_samples = len(counts_per_sample)
        lines.append('Num samples: %s' % str(num_samples))
        lines.append('Num observations: %s' % str(num_observations))
        if not num_observations:
            total_count = sum(counts_per_sample_values)
            lines.append('Total count: %s' % str(total_count))
            lines.append('Table density (fraction of non-zero values): %1.4f' % \
                  table.getTableDensity())
        if not suppress_md5:
            lines.append('Table md5 (unzipped): %s' % safe_md5(biom_open(input_fp,'U')))
        lines.append('')

        if qualitative:
            lines.append('Observations/sample summary:')
        else:
            lines.append('Counts/sample summary:')
        lines.append(' Min: %s' % str(min_counts))
        lines.append(' Max: %s' % str(max_counts))
        lines.append(' Median: %s' % str(median_counts))
        lines.append(' Mean: %s' % str(mean_counts))
        lines.append(' Std. dev.: %s' % (str(std(counts_per_sample_values))))
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
            lines.append(' %s: %s' % (k,str(v)))
        
        result['biom-summary'] = lines
        return result 

    def _get_parameters(self):

        return [Parameter(Name='table',
                          Required=True,
                          Type=Table,
                          Help='the input BIOM table'),
                Parameter(Name='num-observations',
                          Required=False,
                          Type=bool,
                          Help=('Counts are presented as number of unique '
                                'observation ids per sample, rather than '
                                'counts of observations per sample'),
                          Default=False)]

CommandConstructor = PrintBiomTableSummary
