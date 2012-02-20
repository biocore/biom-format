#!/usr/bin/env python

from __future__ import division
from biom.table import DenseOTUTable, SparseOTUTable, table_factory, \
        to_sparsedict
from numpy import zeros
import json

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2012, BIOM-Format Project"
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
                conv_data = zeros(json_table['shape'],dtype=float)
                for r,c,v in json_table['data']:
                    conv_data[r,c] = v
                json_table['data'] = [list(row) for row in conv_data]

            else:
                raise ValueError, "Unknown matrix_type"

    if f is None:
        raise ValueError, 'Unknown table type'

    return f(json_table, constructor)
