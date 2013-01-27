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
from biom.util import compute_counts_per_sample_stats
from biom.parse import parse_biom_table

usage = "usage: Detailed usage examples can be found here: http://biom-format.org/documentation/print_biom_table_summary.html"
desc = "Script to summarize sample or observation data in BIOM-formatted files."

parser = OptionParser(usage=usage, description=desc, version=__version__)
parser.set_defaults(verbose=True)

req_group = OptionGroup(parser, 'Required Options')
req_options = [make_option('-i','--input_fp',help='the input BIOM filepath')]
req_group.add_options(req_options)
parser.add_option_group(req_group)

opt_group = OptionGroup(parser, 'Optional Options')
opt_options = [make_option('-m','--sample_mapping_fp',type="string",
                    help='The sample mapping filepath (will add sample metadata to '+\
                    'biom file, if provided) [default: %default]'),
               make_option('--num_observations',action='store_true',
                   help=('Counts are presented as number of unique observation ids '
                         'per sample, rather than counts of observations per '
                         'sample [default: %default]'),default=False),
               ]
opt_group.add_options(opt_options)
parser.add_option_group(opt_group)


def main():
    opts,args = parser.parse_args()
    input_fp = opts.input_fp
    table = parse_biom_table(open(input_fp,'U'))
    min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
     compute_counts_per_sample_stats(table, opts.num_observations)
    num_otus = len(table.ObservationIds)
    
    counts_per_sample_values = counts_per_sample.values()
    
    try:
        sample_md_keys = table.SampleMetadata[0].keys()
    except TypeError:
        sample_md_keys = ["None provided"]
    try:
        observation_md_keys = table.ObservationMetadata[0].keys()
    except TypeError:
        observation_md_keys = ["None provided"]
    
    num_samples = len(counts_per_sample)
    print 'Num samples: %s' % str(num_samples)
    print 'Num otus: %s' % str(num_otus)
    if not opts.num_observations:
        num_observations = sum(counts_per_sample_values)
        print 'Num observations (sequences): %s' % str(num_observations)
        print 'Table density (fraction of non-zero values): %1.4f' % \
              table.getTableDensity()
    print

    if opts.num_observations:
        print 'OTUs/sample summary:'
    else:
        print 'Seqs/sample summary:' 
    print ' Min: %s' % str(min_counts)
    print ' Max: %s' % str(max_counts)
    print ' Median: %s' % str(median_counts)
    print ' Mean: %s' % str(mean_counts)
    print ' Std. dev.: %s' % (str(std(counts_per_sample_values)))
    print ' Sample Metadata Categories: %s' % '; '.join(sample_md_keys)
    print ' Observation Metadata Categories: %s' % '; '.join(observation_md_keys)
     
    print ''
    if opts.num_observations:
        print 'OTUs/sample detail:'
    else:
        print 'Seqs/sample detail:'
    sorted_counts_per_sample = [(v,k) for k,v in counts_per_sample.items()]
    sorted_counts_per_sample.sort()
    total_count = 0
    for v,k in sorted_counts_per_sample:
        total_count += v
        print ' %s: %s' % (k,str(v))

if __name__ == "__main__":
    main()
