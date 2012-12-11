#!/usr/bin/env python
# File created on 20 Dec 2011
from __future__ import division

from optparse import make_option, OptionParser, OptionGroup
from biom.parse import parse_biom_table, parse_mapping, generatedby

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The BIOM-Format project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

OBS_META_TYPES = {'sc_separated': lambda x: [e.strip() for e in x.split(';')],
                  'naive': lambda x: x
                  }
OBS_META_TYPES['taxonomy'] = OBS_META_TYPES['sc_separated']

usage = "usage: Detailed usage examples can be found here: http://biom-format.org/documentation/metadata_addition.html"
desc = "Script to add sample and/or observation metadata to BIOM-formatted files."

parser = OptionParser(usage=usage, description=desc, version=__version__)
parser.set_defaults(verbose=True)

req_group = OptionGroup(parser, 'Required Options')
req_options = [make_option('-i','--input_fp',help='the input BIOM filepath'),
               make_option('-o','--output_fp',help='the output BIOM filepath')]
req_group.add_options(req_options)
parser.add_option_group(req_group)

opt_group = OptionGroup(parser, 'Optional Options')
opt_options = [make_option('-m','--sample_mapping_fp',type="string",
                    help='The sample mapping filepath (will add sample metadata to '+\
                    'biom file, if provided) [default: %default]'),
                    
               make_option('--observation_mapping_fp',type="string",
                    help='The observation mapping filepath (will add observation metadata '+ \
                            'to biom file, if provided) [default: %default]'),
               ]
opt_group.add_options(opt_options)
parser.add_option_group(opt_group)

def main():
    opts,args = parser.parse_args()

    if opts.input_fp is None:
        parser.print_help()
        parser.error('Must specify an input file!')
    if opts.output_fp is None:
        parser.print_help()
        parser.error('Must specify an output file!')
    
    sample_mapping_fp = opts.sample_mapping_fp
    obs_mapping_fp = opts.observation_mapping_fp
    
    if sample_mapping_fp != None:
        sample_mapping = parse_mapping(open(sample_mapping_fp,'U'))
    else:
        sample_mapping = None
    
    if obs_mapping_fp != None:
        process_fns = {'taxonomy': lambda x: [e.strip() for e in x.split(';')]}
        obs_mapping = parse_mapping(open(obs_mapping_fp, 'U'),process_fns=process_fns)
    else:
        obs_mapping = None
    
    if sample_mapping == None and obs_mapping == None:
        parser.print_help()
        parser.error('Must specify sample_mapping and/or obs_mapping.')
    
    table = parse_biom_table(open(opts.input_fp,'U'))
    output_f = open(opts.output_fp,'w')
    
    if sample_mapping:
        table.addSampleMetadata(sample_mapping)
    
    if obs_mapping:
        table.addObservationMetadata(obs_mapping)

    output_f.write(table.getBiomFormatJsonString(generatedby()))
    
    output_f.close()
    
if __name__ == "__main__":
    main()
