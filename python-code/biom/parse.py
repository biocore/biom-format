#!/usr/bin/env python
# File created on 19 dec 2011
from __future__ import division
from biom.table import DenseOTUTable, SparseOTUTable, table_factory, \
        to_sparsedict
import json
import numpy
from string import strip

# to gut
from qiime.parse import (process_otu_table_sample_ids, parse_otu_table, 
                         parse_mapping_file_to_dict)

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2007-2012, BIOM Format"
__credits__ = ["Justin Kuczynski", "Daniel McDonald", "Greg Caporaso"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "0.9dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Development"

MATRIX_ELEMENT_TYPE = {'int':int,'float':float,'str':str,
                       u'int':int,u'float':float,u'str':str}

def parse_biom_table(json_fh,constructor=None):
    """parses a biom format otu table into a rich otu table object

    input is an open filehandle or compatable object (e.g. list of lines)

    sparse/dense will be determined by "matrix_type" in biom file, and 
    either a SparseOTUTable or DenseOTUTable object will be returned
    note that sparse here refers to the compressed format of [row,col,count]
    dense refers to the full / standard matrix representations
    """
    return parse_biom_table_str(''.join(json_fh),constructor=constructor)

def parse_biom_otu_table(json_table, constructor=None):
    """Parse a biom otu table type

    Constructor must have a _biom_type of "otu table"
    """
    if constructor is None:
        if json_table['matrix_type'].lower() == 'sparse':
            constructor = SparseOTUTable
        elif json_table['matrix_type'].lower() == 'dense':
            constructor = DenseOTUTable
        else:
            raise ValueError, "Unknown matrix_type"
    else:
        if constructor._biom_type.lower() != 'otu table':
            raise ValueError, "constructor must subclass OTUTable"

    sample_ids = [col['id'] for col in json_table['columns']]
    # null metadata -> None object in metadata list 
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]

    table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                              sample_metadata, obs_metadata, 
                              constructor=constructor)

    return table_obj

# map table types -> parsing methods
BIOM_TYPES = {'otu table':parse_biom_otu_table}
def parse_biom_table_str(json_str,constructor=None):
    """Parses a JSON string of the Biom table into a rich table object.
   
    If constructor is none, the constructor is determined based on BIOM
    information
    """
    json_table = json.loads(json_str)

    if constructor is None:
        f = BIOM_TYPES.get(json_table['type'].lower(), None)
    else:
        f = BIOM_TYPES.get(constructor._biom_type.lower(), None)

        # ugh.
        # convert matrix data if the biom type doesn't match matrix type
        # of the table objects
        if constructor._biom_matrix_type != json_table['matrix_type'].lower():
            if json_table['matrix_type'] == 'dense':
                # dense -> sparse
                conv_data = []
                for row_idx,row in enumerate(json_table['data']):
                    for col_idx, value in enumerate(row):
                        if value == 0:
                            continue
                        conv_data.append([row_idx,col_idx,value])
                json_table['data'] = conv_data

            elif json_table['matrix_type'] == 'sparse':
                # sparse -> dense
                conv_data = numpy.zeros(json_table['shape'],dtype=float)
                for r,c,v in json_table['data']:
                    conv_data[r,c] = v
                json_table['data'] = [list(row) for row in conv_data]

            else:
                raise ValueError, "Unknown matrix_type"

    if f is None:
        raise ValueError, 'Unknown table type'

    return f(json_table, constructor)

def parse_otu_table_to_rich_otu_table(lines,count_map_f=int,dense=False,mapping_f=None):
    """parses an otu table (tab delimited) (sample ID x OTU ID map)

    Returns a rich otu table object (a subclass of Table), 
    sparse by default (or see parameter 'dense')
    """
    if mapping_f != None:
        sample_metadata_d = parse_mapping_file_to_dict(mapping_f)[0]
    else:
        sample_metadata_d = None
    
    if dense:
        sample_ids, otu_ids, otu_table, metadata = parse_otu_table(lines,count_map_f=count_map_f)
        if len(metadata) > 0:
            metadata = [{'taxonomy':elem} for elem in metadata]
        else:
            metadata = None
        if sample_metadata_d != None:
            sample_metadata = [sample_metadata_d[sample_id] for sample_id in sample_ids]
        else:
            sample_metadata = None
        table_obj = DenseOTUTable(Data=otu_table,
        SampleIds=sample_ids, ObservationIds=otu_ids,
        SampleMetadata=sample_metadata, ObservationMetadata=metadata)
        return table_obj

    otu_ids = []
    metadata = []
    sample_ids = []
    # iterate over lines in the OTU table -- keep track of line number 
    # to support legacy (Qiime 1.2.0 and earlier) OTU tables
    two_d_dict = {} # {(row,col):value}, row is otu/observaiton
    otu_idx = 0 # keep track of observation/otu lines, (skip comments and headers, etc.)
    for i, line in enumerate(lines):
        line = line.strip()
        if line:
            if i == 1 and line.startswith('#OTU ID') and not sample_ids:
                # we've got a legacy OTU table
                try:
                    sample_ids, has_metadata = process_otu_table_sample_ids(
                     line.strip().split('\t')[1:])
                except ValueError:
                    raise ValueError, \
                     "Error parsing sample IDs in OTU table. Appears to be a"+\
                     " legacy OTU table. Sample ID line:\n %s" % line
            elif not line.startswith('#'):
                if not sample_ids:
                    # current line is the first non-space, non-comment line 
                    # in OTU table, so contains the sample IDs
                    try:
                        sample_ids, has_metadata = process_otu_table_sample_ids(
                         line.strip().split('\t')[1:])
                    except ValueError:
                        raise ValueError,\
                         "Error parsing sample IDs in OTU table."+\
                         " Sample ID line:\n %s" % line
                else:
                    # current line is OTU line in OTU table
                    fields = line.split('\t')

                    # if there is OTU metadata the last column gets appended
                    # to the metadata list
                    # otherwise all columns are appended to otu_table
                    if has_metadata:
                        abund_fields = fields[1:-1]
                        metadata.append({'taxonomy':map(strip, fields[-1].split(';'))})
                    else:
                        abund_fields = fields[1:]
                        metadata.append(None)

                    # added in a try/except to handle OTU tables containing
                    # floating numbers
                    try:
                        abunds = map(count_map_f,abund_fields)
                    except ValueError:
                        abunds = map(float,abund_fields)

                    if all([abund==0 for abund in abunds]):
                        continue #don't increment otu_idx counter
                    else:
                        # grab the OTU ID
                        otu_id = fields[0].strip()
                        otu_ids.append(otu_id)

                    for j, abund in enumerate(abunds):
                        if abund != 0:
                            two_d_dict[(otu_idx,j)] = abund
                    otu_idx += 1 # this is needed for indexing into two_d_dict
                    # this sets dimensions of matrix, so it must be accurate

    if sample_metadata_d != None:
        sample_metadata = [sample_metadata_d[sample_id] for sample_id in sample_ids]
    else:
        sample_metadata = None

    data = to_sparsedict(two_d_dict)
    table_obj = SparseOTUTable(Data=data, 
        SampleIds=sample_ids, ObservationIds=otu_ids,
        SampleMetadata=sample_metadata, ObservationMetadata=metadata)
        
    return(table_obj)

def convert_otu_table_to_biom(otu_table_f,count_map_f=int,dense=False,mapping_f=None):
    """ """
    otu_table = parse_otu_table_to_rich_otu_table(otu_table_f,
                                                  count_map_f=count_map_f,
                                                  dense=dense,
                                                  mapping_f=mapping_f)
    return otu_table.getBiomFormatJsonString()

def convert_biom_to_otu_table(biom_f):
    """ """
    otu_table = parse_biom_table(biom_f)
    if otu_table.ObservationMetadata != None and \
       'taxonomy' in otu_table.ObservationMetadata[0]:
        return otu_table.delimitedSelf(header_key='taxonomy', 
                                       header_value='Consensus Lineage',
                                       metadata_formatter=lambda x: ';'.join(x))
    else:        
        return otu_table.delimitedSelf()
