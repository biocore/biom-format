#!/usr/bin/env python
# File ported from QIIME on 26 Jan 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The BIOM-Format project"
__credits__ = ["Greg Caporaso", "Daniel McDonald","Justin Kuczynski",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from optparse import make_option, OptionParser, OptionGroup
from numpy import std
from biom.util import (compute_counts_per_sample_stats, 
                       safe_md5,
                       biom_open)
from biom.parse import parse_biom_table

usage = "usage: Detailed usage examples can be found here: http://biom-format.org/documentation/summarizing_biom_tables.html"
desc = "Script to summarize sample or observation data in BIOM-formatted files."

parser = OptionParser(usage=usage, description=desc, version=__version__)
parser.set_defaults(verbose=True)

req_group = OptionGroup(parser, 'Required Options')
req_options = [make_option('-i','--input_fp',help='the input BIOM filepath')]
req_group.add_options(req_options)
parser.add_option_group(req_group)

opt_group = OptionGroup(parser, 'Optional Options')
opt_options = [make_option('-o','--output_fp',type="string",default=None,
                    help='Path to write output summary [default: write to stdout]'),
               make_option('--num_observations',action='store_true',
                   help=('Counts are presented as number of unique observation ids '
                         'per sample, rather than counts of observations per '
                         'sample [default: %default]'),default=False),
               make_option('--suppress_md5',action='store_true',
                   help=('Do not include the md5sum of the table in the output' 
                   ' (can be useful if you\'re concerned about runtime of this script)'
                   ' [default: %default]'),default=False),
               ]
opt_group.add_options(opt_options)
parser.add_option_group(opt_group)


def main():
    opts,args = parser.parse_args()

    if opts.input_fp is None:
        parser.print_help()
        parser.error('Must specify an input file!')
        
    input_fp = opts.input_fp
    output_fp = opts.output_fp
    table = parse_biom_table(biom_open(input_fp,'U'))
    min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
     compute_counts_per_sample_stats(table, opts.num_observations)
    num_observations = len(table.ObservationIds)
    suppress_md5 = opts.suppress_md5
    
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
    if not opts.num_observations:
        total_count = sum(counts_per_sample_values)
        lines.append('Total count: %s' % str(total_count))
        lines.append('Table density (fraction of non-zero values): %1.4f' % \
              table.getTableDensity())
    if not suppress_md5:
        lines.append('Table md5 (unzipped): %s' % safe_md5(biom_open(input_fp,'U')))
    lines.append('')

    if opts.num_observations:
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
    if opts.num_observations:
        lines.append('Observations/sample detail:')
    else:
        lines.append('Counts/sample detail:')
    sorted_counts_per_sample = [(v,k) for k,v in counts_per_sample.items()]
    sorted_counts_per_sample.sort()
    for v,k in sorted_counts_per_sample:
        lines.append(' %s: %s' % (k,str(v)))
    
    if output_fp != None:
        open(output_fp,'w').write('\n'.join(lines))
    else:
        print '\n'.join(lines)

if __name__ == "__main__":
    main()
