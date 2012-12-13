#!/usr/bin/env python

"""Validate a Biological Observation Matrix (biom) formatted file

For more details, go to: http://biom-format.org
"""

import json
from httplib import HTTP 
from urlparse import urlparse 
from operator import and_
import dateutil.parser

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jose Clemente", "Greg Caporaso", 
               "Jai Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", 
               "Sue Huse"]
__url__ = "http://biom-format.org"
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

def is_int(x):
    """Return True if x is int"""
    return isinstance(x, int)

def valid_format_url(table):
    """Check if the format_url is correct"""
    if VERBOSE:
        print "Validating format_url..."
        
    if table['format_url'] != FORMAT_URL:
        raise ValueError, "Invalid format_url"

def valid_shape(table):
    """A matrix header is (int, int) representing the size of a 2D matrix"""
    if VERBOSE:
        print "Validating shape..."
    
    a,b = table['shape']

    if not (is_int(a) and is_int(b)):
        raise ValueError, "'shape' values do not appear to be integers"

def valid_matrix_type(table):
    """Returns True if x is a valid matrix type"""
    if VERBOSE:
        print "Validating matrix_type..."
        
    if table['matrix_type'] not in MATRIX_TYPES:
        raise ValueError, "Unknown matrix_type"

def valid_matrix_element_type(table):
    """Return True if table['matrix_element_type'] is a valid element type"""
    if VERBOSE:
        print "Validating matrix_element_type..."
        
    if table['matrix_element_type'] not in ELEMENT_TYPES:
        raise ValueError, "Unknown matrix_element_type"

def valid_datetime(table):
    """Verify datetime can be parsed

    Expects ISO 8601 datetime format (for example, 2011-12-19T19:00:00
                                      note that a 'T' separates the date 
                                      and time)
    """
    if VERBOSE:
        print "Validating datetime..."
        
    try:
        foo = dateutil.parser.parse(table['date'])
    except:
        raise ValueError, "Timestamp does not appear to be ISO 8601"

def valid_sparse_data(table):
    """All index positions must be integers and values are of dtype"""
    if VERBOSE:
        print "Validating data (sparse)..."
        
    dtype = ELEMENT_TYPES[table['matrix_element_type']]
    n_rows, n_cols = table['shape']
    n_rows -= 1 # adjust for 0-based index
    n_cols -= 1 # adjust for 0-based index

    for idx, coord in enumerate(table['data']):
        try:
            x,y,val = coord
        except:
            raise ValueError, "Bad matrix entry idx %d: %s" % (idx,repr(coord))

        if not is_int(x) or not is_int(y):
            raise ValueError, "Bad x or y type at idx %d: %s" % (idx,repr(coord))

        if not isinstance(val, dtype):
            raise ValueError, "Bad value at idx %d: %s" % (idx,repr(coord))

        if x < 0 or x > n_rows:
            raise ValueError, "x out of bounds at idx %d: %s" % (idx,repr(coord))

        if y < 0 or y > n_cols:
            raise ValueError, "y out of bounds at idx %d: %s" % (idx,repr(coord))

def valid_dense_data(table):
    """All elements must be of dtype and correspond to shape"""
    if VERBOSE:
        print "Validating data (dense)..."
        
    dtype = ELEMENT_TYPES[table['matrix_element_type']]
    n_rows, n_cols = table['shape']
   
    for row in table['data']:
        if len(row) != n_cols:
            raise ValueError, "Incorrect number of cols: %s" % repr(row)
        
        if not reduce(and_, [isinstance(v, dtype) for v in row]):
            raise ValueError, "Bad datatype in row: %s" % repr(row)

    if len(table['data']) != n_rows:
        raise ValueError, "Incorrect number of rows in matrix"

def valid_format(table):
    """Format must be the expected version"""
    if VERBOSE:
        print "Validating format..."
        
    if table['format'] != FORMAT_STRING:
        raise ValueError, "Invalid 'format' %s, must be %s" % \
                (table['format'], FORMAT_STRING)

def valid_type(table):
    """Table must be a known table type"""
    if VERBOSE:
        print "Validating type..."
        
    if table['type'].lower() not in BIOM_TYPES:
        raise ValueError, "Unknown BIOM type: %s" % table['type']

def valid_biom(table, check_url=False):
    """Validate a BIOM object
    
    Raises AttributeError if an expected key is missing
    Raises ValueError if the values at a key appear to be malformed
    """
    if VERBOSE:
        print "Validating biom object..."
        
    required_keys = [('format', valid_format),
                     ('format_url', valid_format_url),
                     ('type', valid_type),
                     ('rows', valid_rows),
                     ('columns', valid_columns),
                     ('shape', valid_shape),
                     ('data', valid_data),
                     ('matrix_type', valid_matrix_type),
                     ('matrix_element_type', valid_matrix_element_type),
                     ('generated_by', valid_generated_by),
                     ('id', valid_nullable_id),
                     ('date', valid_datetime)]

    for key,method in required_keys:
        if key not in table:
            raise AttributeError, "MISSING FIELD: '%s'" % key
        method(table)
        #    raise ValueError, "FIELD '%s' INVALID: %s" % (key, repr(table[key]))

    if len(table['rows']) != table['shape'][0]:
        raise ValueError, "Number of rows in 'rows' is not equal to 'shape'"

    if len(table['columns']) != table['shape'][1]:
        raise ValueError, "Number of columns in 'columns' is not equal to 'shape'"

def valid_generated_by(table):
    """Validate the generated_by field"""
    if VERBOSE:
        print "Validating generated_by..."
        
    if not table['generated_by']:
        raise ValueError, "'generated_by' is not populated"
    if not isinstance(table['generated_by'], unicode):
        raise ValueError, "'generated_by' is not a string"

def valid_nullable_id(table):
    """Validate the table id"""
    # this is nullable and don't actually care what is in here
    return

def valid_id(record):
    """Validate id for a row or column"""
    if not record['id']:
        raise ValueError, "'id' in %s appears empty" % record

def valid_metadata(record):
    """Validate the metadata field for a row or column"""
    # this is nullable and don't actually care what is in here
    if record['metadata'] is None:
        return
    if isinstance(record['metadata'], dict):
        return

    raise ValueError, "metadata is neither null or an object"

def valid_rows(table):
    """Validate the 'rows' under 'table'
    
    Raises AttributeError if an expected key is missing
    Raises ValueError if the values at a key appear to be malformed
    """
    if VERBOSE:
        print "Validating rows..."
        
    required_keys = [('id', valid_id), ('metadata', valid_metadata)]
    required_by_type = {}
    required_keys.extend(required_by_type.get(table['type'].lower(), []))

    for idx,row in enumerate(table['rows']):
        for key, method in required_keys:
            if key not in row:
                raise AttributeError, "ROW IDX %d MISSING '%s' FIELD" % (idx,key)
            method(row)

def valid_columns(table):
    """Validate the 'columns' under 'table'

    Raises AttributeError if an expected key is missing
    Raises ValueError if the values at a key appear to be malformed
    """
    if VERBOSE:
        print "Validating columns..."

    required_keys = [('id', valid_id), ('metadata', valid_metadata)]
    required_by_type = {} 
    required_keys.extend(required_by_type.get(table['type'].lower(), []))

    for idx, col in enumerate(table['columns']):
        for key, method in required_keys:
            if key not in col:
                raise AttributeError, "COL IDX %d MISSING '%s' FIELD" % (idx,key)
            method(col)

def valid_data(table):
    """Validate the 'matrix' under 'table'

    Raises AttributeError if an expected key is missing
    Raises ValueError if the values at a key appear to be malformed
    """
    if table['matrix_type'].lower() == 'sparse':
        valid_sparse_data(table)
    elif table['matrix_type'].lower() == 'dense':
        valid_dense_data(table)
    else:
        raise AttributeError, "Unknown matrix type"

try:
    from cogent.util.option_parsing import parse_command_line_parameters, make_option
    cogent_cl_parsing = True
except ImportError:
    from sys import argv
    cogent_cl_parsing = False

FORMAT_URL = "http://biom-format.org"
FORMAT_STRING = "Biological Observation Matrix 1.0.0"
BIOM_TYPES = set(['otu table', 'pathway table', 'function table', 
                  'ortholog table', 'gene table', 'metabolite table', 
                  'taxon table'])
MATRIX_TYPES = set(['sparse', 'dense'])
ELEMENT_TYPES = {'int':int,'str':str,'float':float, 'unicode':unicode}
VERBOSE = False

if cogent_cl_parsing:
    script_info = {}
    script_info['brief_description'] = "Test a biom file for adherence to the format specification."
    script_info['script_description'] = "Test a biom file for adherence to the format specification. This specification is defined at http://biom-format.org."
    script_info['script_usage'] = [("","Validate the my_data.biom file.","%prog -i my_data.biom")]
    script_info['output_description']= ""
    script_info['required_options'] = [
     make_option('-i','--biom_fp',type="existing_filepath",
                 help='the BIological Observation Matrix filepath to validate'),
    ]
    script_info['optional_options'] = [
     make_option('-f','--format-version',type="string", 
             default=FORMAT_STRING,
             help='The specific format string, defaults to [default: %default]')]
    script_info['version'] = __version__

if __name__ == '__main__':
    if cogent_cl_parsing:
        option_parser, opts, args =\
         parse_command_line_parameters(**script_info)
        biom_fp = opts.biom_fp
        VERBOSE = opts.verbose
        FORMAT_STRING = opts.format_version
        valid_biom(json.load(open(biom_fp)))
    else:
        if '-v' in argv:
            VERBOSE = True
            argv.remove('-v')
        elif '--verbose' in argv:
            VERBOSE = True
            argv.remove('--verbose')
            
        if len(argv) == 3:
            biom_fp = open(argv[2])
        elif len(argv) == 5:
            biom_fp = open(argv[2])
            FORMAT_STRING = argv[4]
        else:
            print "Error parsing command.\nUSAGE: biom_validator.py -i my_file.biom"
            print "\nOptional:\n\t-f\tSpecify the format string, default to '%s'" % FORMAT_STRING
            exit()
        valid_biom(json.load(biom_fp))

