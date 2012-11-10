#!/usr/bin/env python

from __future__ import division
from biom.exception import BiomParseException
from biom.table import SparseOTUTable, DenseOTUTable, SparsePathwayTable, \
        DensePathwayTable, SparseFunctionTable, DenseFunctionTable, \
        SparseOrthologTable, DenseOrthologTable, SparseGeneTable, \
        DenseGeneTable, SparseMetaboliteTable, DenseMetaboliteTable,\
        SparseTaxonTable, DenseTaxonTable, table_factory, to_sparse,\
        nparray_to_sparseobj, SparseObj
from biom.csmat import CSMat
import json
from numpy import zeros, asarray, uint32, float64
from string import strip

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Justin Kuczynski", "Daniel McDonald", "Greg Caporaso", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.0.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

MATRIX_ELEMENT_TYPE = {'int':int,'float':float,'unicode':unicode,
                      u'int':int,u'float':float,u'unicode':unicode}

QUOTE = '"'
JSON_OPEN = set(["[", "{"])
JSON_CLOSE = set(["]", "}"])
JSON_SKIP = set([" ","\t","\n",","])
JSON_START = set(["0","1","2","3","4","5","6","7","8","9","{","[",'"'])
def direct_parse_key(biom_str, key):
    """Returns key:value from the biom string, or ""
    
    This method pulls an arbitrary key/value pair out from a BIOM string
    """
    base_idx = biom_str.find('"%s":' % key)
    if base_idx == -1:
        return ""
    else:
        start_idx = base_idx + len(key) + 3 # shift over "key":

    # find the start token
    cur_idx = start_idx 
    while biom_str[cur_idx] not in JSON_START:  
        cur_idx += 1
    
    if biom_str[cur_idx] not in JSON_OPEN:
        # do we have a number?
        while biom_str[cur_idx] not in [",","{","}"]:
            cur_idx += 1
            
    else:
        # we have an object
        stack = [biom_str[cur_idx]]
        cur_idx += 1
        while stack:
            cur_char = biom_str[cur_idx]
            
            if cur_char == QUOTE:
                if stack[-1] == QUOTE:
                    stack.pop()
                else:
                    stack.append(cur_char)
            elif cur_char in JSON_CLOSE:
                try:
                    stack.pop()
                except IndexError: # got an int or float?
                    cur_idx -= 1
                    break
            elif cur_char in JSON_OPEN:
                stack.append(cur_char)
            cur_idx += 1

    return biom_str[base_idx:cur_idx]

def direct_slice_data(biom_str, to_keep, axis):
    """Pull out specific slices from a BIOM string

    biom_str : JSON-formatted BIOM string
    to_keep  : indices to keep
    axis     : either 'samples' or 'observations'
    
    Will raise IndexError if the inices are out of bounds. Fully zerod rows
    or columns are possible and this is _not_ checked. 
    """
    if axis not in ['observations','samples']:
        raise IndexError, "Unknown axis type"

    # it would be nice if all of these lookups could be done in a single
    # traversal of biom_str, but it likely is at the cost of code complexity
    shape_kv_pair = direct_parse_key(biom_str, "shape")
    if shape_kv_pair == "":
        raise ValueError, "biom_str does not appear to be in BIOM format!"

    data_fields = direct_parse_key(biom_str, "data")
    if data_fields == "":
        raise ValueError, "biom_str does not appear to be in BIOM format!"

    matrix_type_kv_pair = direct_parse_key(biom_str, "matrix_type")
    if matrix_type_kv_pair == "":
        raise ValueError, "biom_str does not appear to be in BIOM format!"

    # determine shape
    raw_shape = shape_kv_pair.split(':')[-1].replace("[","").replace("]","")
    n_rows, n_cols = map(int, raw_shape.split(","))
   
    # slice to just data
    data_start = data_fields.find('[') + 1 # find [[ for start, and move to [
    data_fields = data_fields[data_start:len(data_fields)-1] # trim trailing ]
    
    # determine matrix type
    matrix_type = matrix_type_kv_pair.split(':')[-1].strip()

    # bounds check
    if min(to_keep) < 0:
        raise IndexError, "Observations to keep are out of bounds!"
    
    # more bounds check and set new shape
    new_shape = "[%d, %d]"
    if axis == 'observations':
        if max(to_keep) >= n_rows:
            raise IndexError, "Observations to keep are out of bounds!"
        new_shape = new_shape % (len(to_keep), n_cols)
    elif axis == 'samples':
        if max(to_keep) >= n_cols:
            raise IndexError, "Samples to keep are out of bounds!"
        new_shape = new_shape % (n_rows, len(to_keep))

    to_keep = set(to_keep)
    new_data = []
    ### decompose these???
    if matrix_type == '"dense"':
        if axis == 'observations':
            # iterate rows, keep those specified
            for row_count, row in enumerate(data_fields.split('],')):
                if row_count in to_keep:
                    new_data.append(row.strip())
        elif axis == 'samples':
            # iterate rows
            for row in data_fields.split('],'):
                new_row = []
                # dive into the cols and keep those specified
                for col_idx,v in enumerate(row.split(',')):
                    if col_idx in to_keep:
                        v = v.strip()
                        new_row.append(v)
                new_data.append("[%s" % ','.join(new_row))
    elif matrix_type == '"sparse"':
        if axis == 'observations':
            # interogate all the datas
            for rcv in data_fields.split('],'):
                r,c,v = rcv.split(',')
                r = r.strip()[1:] # trim off [
                if int(r) in to_keep:
                    new_data.append(rcv.strip())
        elif axis == 'samples':
            # could do sparse obs/samp in one forloop, but then theres the 
            # expense of the additional if-statement in the loop
            for rcv in data_fields.split('],'):
                r,c,v = rcv.split(',')
                if int(c) in to_keep:
                    new_data.append(rcv.strip())
    else:
        print matrix_type, matrix_type == '"dense"'
        raise ValueError, "biom_str does not appear to be in BIOM format!"

    return '"data": [%s]], "shape": %s' % ('],'.join(new_data), new_shape)

def get_axis_indices(biom_str, to_keep, axis):
    """Returns the indices for the associated ids to keep

    biom_str : a BIOM formatted JSON string
    to_keep  : a list of IDs to get indices for
    axis     : either 'samples' or 'observations'
    
    Raises KeyError if unknown key is specified
    """
    to_keep = set(to_keep)
    if axis == 'observations':
        axis_key = 'rows'
        axis_data = direct_parse_key(biom_str, axis_key)
    elif axis == "samples":
        axis_key = 'columns'
        axis_data = direct_parse_key(biom_str, axis_key)
    else:
        raise ValueError, "Unknown axis!"

    if axis_data == "":
        raise ValueError, "biom_str does not appear to be in BIOM format!"

    axis_data = json.loads("{%s}" % axis_data)

    all_ids = set([v['id'] for v in axis_data[axis_key]])
    if not to_keep.issubset(all_ids):
        raise KeyError, "Not all of the to_keep ids are in biom_str!"

    return [i for i,v in enumerate(axis_data[axis_key]) if v['id'] in to_keep]

def light_parse_biom_csmat(biom_str, constructor):
    """Light-weight BIOM parser for CSMat-Sparse objects

    Constructor must match the loaded table type
    """
    if constructor._biom_matrix_type != "sparse":
        raise AttributeError, "Constructor must be sparse!"

    if not isinstance(SparseObj(1,1), CSMat):
        raise AttributeError, "SparseObj is not CSMat!"
    
    # is data: separated by a space?
    data_start = biom_str.find('"data":')
    if biom_str[data_start + 7] == " ":
        start_idx = data_start + 8
    else:
        start_idx = data_start + 7

    #start_idx = biom_str.find('"data":') + 7 
    #print biom_str[start_idx-10:start_idx+10]
    end_idx = biom_str[start_idx:].find(']]') + start_idx
    data = biom_str[start_idx:end_idx]
    new_s = biom_str[:start_idx]
    new_s += '[[0, 0, 1]]'
    new_s += biom_str[(end_idx + 2):]
    
    # get shape
    start_idx = biom_str.find('"shape":') + 10
    end_idx = biom_str[start_idx:start_idx + 30].find('],') + start_idx
    row, col = map(int, biom_str[start_idx:end_idx].replace('[','').split(', '))
    data_mat = CSMat(row,col)
    
    # if we knew in advance the number of nnz, we could preallocate
    # for storage
    for rec in data.replace('[','').split('],'):
        try:
            row,col,val = rec.split(',')
        except:
            raise TypeError, "Data do not appear sparse!"

        data_mat._coo_rows.append(uint32(row))
        data_mat._coo_cols.append(uint32(col))
        data_mat._coo_values.append(float64(val))
    
    t = parse_biom_table_str(new_s, constructor, data_pump=data_mat)

    return t

def parse_biom_table(json_fh,constructor=None):
    """parses a biom format otu table into a rich otu table object

    input is an open filehandle or compatable object (e.g. list of lines)

    sparse/dense will be determined by "matrix_type" in biom file, and 
    either a SparseOTUTable or DenseOTUTable object will be returned
    note that sparse here refers to the compressed format of [row,col,count]
    dense refers to the full / standard matrix representations
    """
    return parse_biom_table_str(''.join(json_fh),constructor=constructor)

def pick_constructor(mat_type, table_type, constructor, valid_constructors):
    """Make sure constructor is sane, attempt to pick one if not specified

    Excepts valid_constructors to be a list in the order of 
    [SparseTable, DenseTable] in which the objects present must subclass the
    objects respectively (eg [SparseOTUTable, DenseOTUTable])

    We do not require the matrix type to be the same as the constructor if the 
    passed in constructor is not None. The motivation is that there are use
    cases for taking a table stored as dense but loaded as sparse.

    Will raise BiomParseError if input_mat_type appears wrong or if the 
    specified constructor appears to be incorrect
    """
    if constructor is None:
        if mat_type.lower() == 'sparse':
            constructor = valid_constructors[0]
        elif mat_type.lower() == 'dense':
            constructor = valid_constructors[1]
        else:
            raise BiomParseException, "Unknown matrix_type"

    if constructor._biom_type.lower() != table_type.lower():
        raise BiomParseException, "constructor must be a biom %s" % table_type

    return constructor

def parse_biom_otu_table(json_table, constructor=None, data_pump=None):
    """Parse a biom otu table type

    Constructor must have a _biom_type of "otu table"
    """
    table_type = 'otu table'
    mat_type = json_table['matrix_type']
    constructors = [SparseOTUTable, DenseOTUTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor, 
                                  shape=json_table['shape'], 
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor, 
                                  shape=json_table['shape'], 
                                  dtype=dtype)

    return table_obj

def parse_biom_pathway_table(json_table, constructor=None, data_pump=None):
    """Parse a biom pathway table
    
    Constructor must have a _biom_type of "pathway table"
    """
    mat_type = json_table['matrix_type']
    table_type = 'pathway table'
    constructors = [SparsePathwayTable, DensePathwayTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:    
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)

    return table_obj

def parse_biom_function_table(json_table, constructor=None, data_pump=None):
    """Parse a biom function table
    
    Constructor must have a _biom_type of "function table"
    """
    mat_type = json_table['matrix_type']
    table_type = 'function table'
    constructors = [SparseFunctionTable, DenseFunctionTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)

    return table_obj

def parse_biom_ortholog_table(json_table, constructor=None, data_pump=None):
    """Parse a biom ortholog table

    Constructor must have a _biom_type of "ortholog table"
    """
    mat_type = json_table['matrix_type']
    table_type = 'ortholog table'
    constructors = [SparseOrthologTable, DenseOrthologTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)

    return table_obj

def parse_biom_gene_table(json_table, constructor=None, data_pump=None):
    """Parse a biom gene table
    
    Constructor must have a _biom_type of "gene table"
    """
    mat_type = json_table['matrix_type']
    table_type = 'gene table'
    constructors = [SparseGeneTable, DenseGeneTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    return table_obj

def parse_biom_metabolite_table(json_table, constructor=None, data_pump=None):
    """Parse a biom metabolite table

    Constructor must have a _biom_type of "metabolite table"
    """
    mat_type = json_table['matrix_type']
    table_type = 'metabolite table'
    constructors = [SparseMetaboliteTable, DenseMetaboliteTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    return table_obj

def parse_biom_taxon_table(json_table, constructor=None, data_pump=None):
    """Parse a biom taxon table

    Constructor must have a _biom_type of "taxon table"
    """
    mat_type = json_table['matrix_type']
    table_type = 'taxon table'
    constructors = [SparseTaxonTable, DenseTaxonTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor,
                                  shape=json_table['shape'],
                                  dtype=dtype)

    return table_obj

# map table types -> parsing methods
BIOM_TYPES = {'otu table':parse_biom_otu_table,
              'pathway table':parse_biom_pathway_table,
              'function table':parse_biom_function_table,
              'ortholog table':parse_biom_ortholog_table,
              'gene table':parse_biom_gene_table,
              'metabolite table':parse_biom_metabolite_table,
              'taxon table':parse_biom_taxon_table}

def parse_biom_table_str(json_str,constructor=None, data_pump=None):
    """Parses a JSON string of the Biom table into a rich table object.
   
    If constructor is none, the constructor is determined based on BIOM
    information

    data_pump is to allow the injection of a pre-parsed data object
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
                raise BiomParseException, "Unknown matrix_type"

    if f is None:
        raise BiomParseException, 'Unknown table type'

    return f(json_table, constructor, data_pump)

OBS_META_TYPES = {'sc_separated': lambda x: [e.strip() for e in x.split(';')],
                  'naive': lambda x: x
                  }
OBS_META_TYPES['taxonomy'] = OBS_META_TYPES['sc_separated']

def parse_classic_table_to_rich_table(lines, sample_mapping, obs_mapping, process_func,
        constructor, **kwargs):
    """Parses an table (tab delimited) (observation x sample)

    sample_mapping : can be None or {'sample_id':something}
    obs_mapping : can be none or {'observation_id':something}
    """
    sample_ids, obs_ids, data, t_md, t_md_name = parse_classic_table(lines, 
                                                        **kwargs)

    # if we have it, keep it
    if t_md is None:
        obs_metadata = None
    else:
        obs_metadata = [{t_md_name:process_func(v)} for v in t_md]

    if sample_mapping is None:
        sample_metadata = None
    else:
        sample_metadata = [sample_mapping[sample_id] for sample_id in sample_ids]

    # will override any metadata from parsed table
    if obs_mapping is not None:
        obs_metadata = [obs_mapping[obs_id] for obs_id in obs_ids]

    if constructor._biom_matrix_type == 'sparse':
        data = nparray_to_sparseobj(data)
    
    return table_factory(data, sample_ids, obs_ids, sample_metadata, 
                         obs_metadata, constructor=constructor)

def parse_classic_table(lines, delim='\t', dtype=float, header_mark=None, \
        md_parse=None):
    """Parse a classic table into (sample_ids, obs_ids, data, metadata, md_name)

    If the last column does not appear to be numeric, interpret it as 
    observation metadata, otherwise None.

    md_name is the column name for the last column if non-numeric

    NOTE: this is intended to be close to how QIIME classic OTU tables are
    parsed with the exception of the additional md_name field
    """
    if not isinstance(lines, list):
        try:
            lines = lines.readlines()
        except AttributeError:
            raise BiomParseException, "Input needs to support readlines or be indexable"

    # find header, the first line that is not empty and does not start with a #
    for idx,l in enumerate(lines):
        if not l.strip():
            continue
        if not l.startswith('#'):
            break
        if header_mark and l.startswith(header_mark):
            break

    if idx == 0:
        data_start = 1
        header = lines[0].strip().split(delim)[1:]
    else:
        if header_mark is not None:
            data_start = idx + 1
            header = lines[idx].strip().split(delim)[1:]
        else:
            data_start = idx
            header = lines[idx-1].strip().split(delim)[1:]

    # attempt to determine if the last column is non-numeric, ie, metadata
    first_values = lines[data_start].strip().split(delim)
    last_value = first_values[-1]
    last_column_is_numeric = True

    if '.' in last_value:
        try:
            float(last_value)
        except ValueError:
            last_column_is_numeric = False
    else:
        try:
            int(last_value)
        except ValueError:
            last_column_is_numeric = False
    
    # determine sample ids
    if last_column_is_numeric:
        md_name = None
        metadata = None
        samp_ids = header[:]
    else:
        md_name = header[-1]
        metadata = []
        samp_ids = header[:-1]
    
    data = []
    obs_ids = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            continue

        fields = line.strip().split(delim)
        obs_ids.append(fields[0])

        if last_column_is_numeric:
            values = map(dtype, fields[1:])
        else:
            values = map(dtype, fields[1:-1])

            if md_parse is not None:
                metadata.append(md_parse(fields[-1]))
            else:
                metadata.append(fields[-1])

        data.append(values)

    return samp_ids, obs_ids, asarray(data), metadata, md_name

def parse_mapping(lines, strip_quotes=True, suppress_stripping=False):
    """Parser for map file that relates samples to metadata.
    
    Format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    Result: {first_column:{column_i:value}}, where i > 0

    Assumes the first column in the mapping file is the id

    NOTE: code pulled and modified from QIIME (http://qiime.org)
    """
    if hasattr(lines,"upper"):
        # Try opening if a string was passed
        try:
            lines = open(lines,'U')
        except IOError:
            raise BiomParseException,\
             ("A string was passed that doesn't refer "
              "to an accessible filepath.")
    
    if strip_quotes:
        if suppress_stripping:
            # remove quotes but not spaces
            strip_f = lambda x: x.replace('"','')
        else:
            # remove quotes and spaces
            strip_f = lambda x: x.replace('"','').strip()
    else:
        if suppress_stripping:
            # don't remove quotes or spaces
            strip_f = lambda x: x
        else:
            # remove spaces but not quotes
            strip_f = lambda x: x.strip()

    # Create lists to store the results
    mapping_data = []
    header = []
    comments = []

    # Begin iterating over lines
    for line in lines:
        line = strip_f(line)
        if not line or (suppress_stripping and not line.strip()):
            # skip blank lines when not stripping lines
            continue

        if line.startswith('#'):
            line = line[1:]
            if not header:
                header = line.strip().split('\t')
            else:
                comments.append(line)
        else:
            # Will add empty string to empty fields
            tmp_line = map(strip_f, line.split('\t'))
            if len(tmp_line)<len(header):
                tmp_line.extend(['']*(len(header)-len(tmp_line)))
            mapping_data.append(tmp_line)

    if not header:
        raise BiomParseException, "No header line was found in mapping file."
    if not mapping_data:
        raise BiomParseException, "No data found in mapping file."

    first_col = [i[0] for i in mapping_data]
    if len(first_col) != len(set(first_col)):
        raise BiomParseException, "First column values are not unique! Cannot be ids."

    mapping = {}
    for vals in mapping_data:
        mapping[vals[0]] = dict([(k,v) for k,v in zip(header[1:], vals[1:])])

    return mapping

def generatedby():
    """Returns a generated by string"""
    return 'BIOM-Format %s' % __version__

def convert_table_to_biom(table_f,sample_mapping, obs_mapping, process_func, constructor,
                          **kwargs):
    """Convert a contigency table to a biom table
    
    sample_mapping : dict of {'sample_id':metadata} or None
    obs_mapping : dict of {'obs_id':metadata} or None
    process_func: a function to transform observation metadata
    constructor : a biom table type
    dtype : type of table data
    """
    otu_table = parse_classic_table_to_rich_table(table_f, sample_mapping, 
                                                  obs_mapping, process_func,
                                                  constructor, **kwargs)
    return otu_table.getBiomFormatJsonString(generatedby())

def convert_biom_to_table(biom_f, header_key=None, header_value=None, \
        md_format=None):
    """Convert a biom table to a contigency table"""
    table = parse_biom_table(biom_f)

    if md_format is None:
        md_format = lambda x: '; '.join(x)

    if table.ObservationMetadata is None:
        return table.delimitedSelf()
    
    if header_key in table.ObservationMetadata[0]:
        return table.delimitedSelf(header_key=header_key, 
                                       header_value=header_value,
                                       metadata_formatter=md_format)
    else:
        return table.delimitedSelf()
