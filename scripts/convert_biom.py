#!/usr/bin/env python
# File created on 20 Dec 2011
from __future__ import division

from optparse import make_option, OptionParser, OptionGroup
from biom.table import SparseOTUTable, DenseOTUTable, SparsePathwayTable, \
        DensePathwayTable, SparseFunctionTable, DenseFunctionTable, \
        SparseOrthologTable, DenseOrthologTable, SparseGeneTable, \
        DenseGeneTable, SparseMetaboliteTable, DenseMetaboliteTable,\
        SparseTaxonTable, DenseTaxonTable, table_factory
from biom.parse import parse_biom_table, parse_mapping, convert_biom_to_table, \
        convert_table_to_biom, generatedby

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The BIOM-Format project"
__credits__ = ["Greg Caporaso", "Daniel McDonald",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

BIOM_TYPES = {'otu table':[SparseOTUTable, DenseOTUTable], 
              'pathway table':[SparsePathwayTable, DensePathwayTable], 
              'function table':[SparseFunctionTable, DenseFunctionTable],
              'ortholog table':[SparseOrthologTable, DenseOrthologTable],
              'gene table':[SparseGeneTable, DenseGeneTable],
              'metabolite table':[SparseMetaboliteTable, DenseMetaboliteTable],
              'taxon table':[SparseTaxonTable, DenseTaxonTable]}

OBS_META_TYPES = {'sc_separated': lambda x: [e.strip() for e in x.split(';')],
                  'naive': lambda x: x
                  }
OBS_META_TYPES['taxonomy'] = OBS_META_TYPES['sc_separated']

usage = "usage: Detailed usage examples can be found here: http://biom-format.org/documentation/biom_conversion.html"
desc = "Script to convert biom formatted files."

parser = OptionParser(usage=usage, description=desc, version=__version__)
parser.set_defaults(verbose=True)

req_group = OptionGroup(parser, 'Required Options')
req_options = [make_option('-i','--input_fp',type="string",
                           help='the input filepath'),
               make_option('-o','--output_fp',type="string",
                           help='the output filepath')]
req_group.add_options(req_options)
parser.add_option_group(req_group)

opt_group = OptionGroup(parser, 'Optional Options')
opt_options = [make_option('-t','--biom_type',type='choice',
                    choices=['sparse','dense'],default='sparse',
                    help="Type of biom file to write (dense or sparse) when "
                          "passed a classic table [default: %default]"),
               make_option('-b','--biom_to_classic_table',
                    action='store_true', help="Convert biom file to classic "
                    "table file [default: convert "
                    "classic table file to biom file]",default=False),
               make_option('--sparse_biom_to_dense_biom',action='store_true',
                    help="Convert sparse biom file to a dense biom file " 
                         "[default: convert "
                         "classic table file to biom file]",default=False),
               make_option('--dense_biom_to_sparse_biom',action='store_true',
                    help="Convert dense biom file to a sparse biom file "
                    "[default: convert "
                    "classic table file to biom file]",default=False),
               make_option('-m','--sample_mapping_fp',type="string",
                    help='The mapping filepath (will add sample metadata to '
                    'biom file, if provided) [default: %default]'),
               make_option('--observation_mapping_fp',type="string",
                    help='The mapping filepath (will add observation metadata '
                    'to biom file, if provided) [default: %default]'),
               make_option('--header_key',type="string",default=None,
                    help='Pull this key from observation metadata within a '
                    'biom file when writing a classic table [default: no '
                    'observation metadata will be written]'),
               make_option('--output_metadata_id',type="string",default=None,
                    help='Name for the observation metadata column when '
                    'writing a biom file as a classic table [default: same '
                    'name as in biom-formatted table]'),
               make_option('--process_obs_metadata',type="choice",
                           choices=OBS_META_TYPES.keys(), default='naive',
                    help='Process metadata associated with observations when '
                    'converting from a classic table. Must be one of: ' +
                    ', '.join(OBS_META_TYPES.keys()) + ' [default: %default]'),
               make_option('--biom_table_type',type="string",default=None,
                    help='The biom table type to get converted into. Required '
                    'when converting a classic table file to biom file. Must '
                    'be one of: %s' % ', '.join(BIOM_TYPES.keys()))]
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

    biom_to_classic_table = opts.biom_to_classic_table
    sparse_biom_to_dense_biom = opts.sparse_biom_to_dense_biom
    dense_biom_to_sparse_biom = opts.dense_biom_to_sparse_biom
    process_obs_metadata = opts.process_obs_metadata
    
    if sum([biom_to_classic_table,
            sparse_biom_to_dense_biom,
            dense_biom_to_sparse_biom]) > 1:
        parser.print_help()
        option_parser.error("The --biom_to_classic_table, --sparse_biom_to_dense_biom, "
         "and --dense_biom_to_sparse_biom options are mutually exclusive. Pass only one at a time.")
    
    input_f = open(opts.input_fp,'U')
    output_f = open(opts.output_fp,'w')
    
    #dense = opts.biom_type == 'dense'
    count_map_f = int
    sample_mapping_fp = opts.sample_mapping_fp
    obs_mapping_fp = opts.observation_mapping_fp
    
    if sample_mapping_fp != None:
        sample_mapping = parse_mapping(open(sample_mapping_fp,'U'))
    else:
        sample_mapping = None
    
    if obs_mapping_fp != None:
        obs_mapping = parse_mapping(open(obs_mapping_fp, 'U'))
    else:
        obs_mapping = None
    
    # if the user does not specify a name for the output metadata column, set it to the
    # same as the header key
    header_key = opts.header_key
    output_metadata_id = opts.output_metadata_id or header_key
    if biom_to_classic_table:
        try:
            output_f.write(convert_biom_to_table(\
                           input_f, header_key, output_metadata_id))
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"
    elif sparse_biom_to_dense_biom:
        try:
            table = parse_biom_table(input_f)
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"        

        conv_constructor = BIOM_TYPES[table._biom_type.lower()][1]
        conv_table = table_factory(table._data, table.SampleIds, 
                        table.ObservationIds, table.SampleMetadata, 
                        table.ObservationMetadata, table.TableId, 
                        constructor=conv_constructor)
        output_f.write(conv_table.getBiomFormatJsonString(generatedby()))
    elif dense_biom_to_sparse_biom:
        try:
            table = parse_biom_table(input_f)
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"

        conv_constructor = BIOM_TYPES[table._biom_type.lower()][0]
        conv_table = table_factory(table._data, table.SampleIds, 
                        table.ObservationIds, table.SampleMetadata, 
                        table.ObservationMetadata, table.TableId, 
                        constructor=conv_constructor)

        output_f.write(conv_table.getBiomFormatJsonString(generatedby()))
    else:
        if opts.biom_table_type is None:
            parser.error('Must specify the BIOM table type: %s' % \
                    ', '.join(BIOM_TYPES.keys()))
        else:
            biom_table_type = opts.biom_table_type.lower()
        if biom_table_type not in BIOM_TYPES:
            parser.error('Unknown BIOM table type, must be one of: %s' % \
                    ', '.join(BIOM_TYPES.keys()))
        if opts.biom_type is None or opts.biom_type not in ['dense', 'sparse']:
            parser.error('Must specify the BIOM matrix type, ' + \
                    'either "dense" or "sparse"')

        idx = 0 if opts.biom_type == 'sparse' else 1
        constructor = BIOM_TYPES[biom_table_type][idx]
        
        try:
            output_f.write(convert_table_to_biom(input_f,sample_mapping, 
                                obs_mapping, OBS_META_TYPES[process_obs_metadata], constructor))
        except ValueError:
            raise ValueError, "Input does not look like a classic table. Do you need to pass -b?"
    input_f.close()
    output_f.close()
    
if __name__ == "__main__":
    main()
