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
               make_option('--hierarchical_fields',type="string",
                    help=('comma-separated list of the metadata '
                          'fields to split on semi-colons. this is useful '
                          'for hierarchical data such as taxonomy or functional '
                          'category [default: %default]'),
                          default=None),
               make_option('--int_fields',type="string",
                    help=('comma-separated list of the metadata '
                          'fields to cast to integers. this is useful '
                          'for integer data such as "DaysSinceStart" [default: %default]'),
                          default=None),
               make_option('--float_fields',type="string",
                    help=('comma-separated list of the metadata '
                          'fields to cast to floating point numbers. this is useful '
                          'for real number data such as "pH" [default: %default]'),
                          default=None),
               make_option('--observation_header',type="string",
                    help=('comma-separated list of the observation metadata '
                          'field names. this is useful if a header line is not '
                          'provided with the metadata, if you want to rename the '
                          'fields, or if you want to include only the first n '
                          'fields where n is the number of entries provided here '
                          ' [default: use header from observation_mapping_fp]'),
                          default=None),
               make_option('--sample_header',type="string",
                    help=('comma-separated list of the sample metadata '
                          'field names. this is useful if a header line is not '
                          'provided with the metadata, if you want to rename the '
                          'fields, or if you want to include only the first n '
                          'fields where n is the number of entries provided here '
                          ' [default: use header from sample_mapping_fp]'),
                          default=None),
               ]
opt_group.add_options(opt_options)
parser.add_option_group(opt_group)

def split_on_semicolons(x):
    return [e.strip() for e in x.split(';')]

def int_(x):
    try:
        return int(x)
    except ValueError:
        return x

def float_(x):
    try:
        return float(x)
    except ValueError:
        return x

def main():
    opts,args = parser.parse_args()

    if opts.input_fp is None:
        parser.print_help()
        parser.error('Must specify an input file!')
    if opts.output_fp is None:
        parser.print_help()
        parser.error('Must specify an output file!')
    
    ## process header information, if provided
    observation_header = opts.observation_header
    sample_header = opts.sample_header
    if opts.observation_header != None:
        observation_header = observation_header.split(',')
    if opts.sample_header != None:
        sample_header = sample_header.split(',')
    
    ## define metadata processing functions, if any
    process_fns = {}
    hierarchical_fields = opts.hierarchical_fields
    if hierarchical_fields != None:
        process_fns.update({}.fromkeys(hierarchical_fields.split(','),
                                      split_on_semicolons))
    
    int_fields = opts.int_fields
    if int_fields != None:
        process_fns.update({}.fromkeys(int_fields.split(','),
                                      int_))
    
    float_fields = opts.float_fields
    if float_fields != None:
        process_fns.update({}.fromkeys(float_fields.split(','),
                                      float_))

    ## parse mapping files
    sample_mapping_fp = opts.sample_mapping_fp
    obs_mapping_fp = opts.observation_mapping_fp
    if sample_mapping_fp != None:
        sample_mapping = parse_mapping(open(sample_mapping_fp,'U'),
                                       process_fns=process_fns,
                                       header=sample_header)
    else:
        sample_mapping = None
    
    if obs_mapping_fp != None:
        obs_mapping = parse_mapping(open(obs_mapping_fp, 'U'),
                                    process_fns=process_fns,
                                    header=observation_header)
    else:
        obs_mapping = None
    
    if sample_mapping == None and obs_mapping == None:
        parser.print_help()
        parser.error('Must specify sample_mapping and/or obs_mapping.')
    
    ## parse the table and open the output file for writing
    output_f = open(opts.output_fp,'w')
    table = parse_biom_table(open(opts.input_fp,'U'))
    
    ## add metadata as necessary
    if sample_mapping:
        table.addSampleMetadata(sample_mapping)
    
    if obs_mapping:
        table.addObservationMetadata(obs_mapping)

    ## write the output file and close it
    output_f.write(table.getBiomFormatJsonString(generatedby()))
    output_f.close()
    
if __name__ == "__main__":
    main()
