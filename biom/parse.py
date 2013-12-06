#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from biom import __version__
from biom.exception import BiomParseException
from biom.table import Table, table_factory, to_sparse,\
        nparray_to_sparseobj, SparseObj
import json
from numpy import zeros, asarray, uint32, float64
from string import strip
import h5py
from scipy.sparse import coo_matrix

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Justin Kuczynski", "Daniel McDonald", "Greg Caporaso", "Jose Carlos Clemente Litran","Morgan Langille"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

MATRIX_ELEMENT_TYPE = {'int':int,'float':float,'unicode':unicode,
                       u'int':int,u'float':float,u'unicode':unicode}

def parse_biom_table(fp):
    if hasattr(fp, 'read'):
        return parse_biom_table_json(json.load(fp))
    else:
        return parse_biom_table_json(json.loads(fp))

def parse_biom_table_json(json_table, data_pump=None):
    """Parse a biom otu table type"""
    constructor = Table
    table_type = 'otu table'
    mat_type = json_table['matrix_type']

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  shape=json_table['shape'], 
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  shape=json_table['shape'], 
                                  dtype=dtype)

    return table_obj

def sc_pipe_separated(x):
    complex_metadata=[]
    for y in x.split('|'):
        simple_metadata=[]
        for e in y.split(';'):
            simple_metadata.append(e.strip())
        complex_metadata.append(simple_metadata)
    return complex_metadata


OBS_META_TYPES = {'sc_separated': lambda x: [e.strip() for e in x.split(';')],
                  'naive': lambda x: x, 'sc_pipe_separated': sc_pipe_separated
                  }
OBS_META_TYPES['taxonomy'] = OBS_META_TYPES['sc_separated']

def parse_classic_table_to_rich_table(lines, sample_mapping, obs_mapping, 
        process_func, **kwargs):
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

    data = nparray_to_sparseobj(data)
    
    return table_factory(data, sample_ids, obs_ids, sample_metadata, 
                         obs_metadata)

def parse_classic_table(lines, delim='\t', dtype=float, header_mark=None, \
        md_parse=None):
    """Parse a classic table into (sample_ids, obs_ids, data, metadata, md_name)

    If the last column does not appear to be numeric, interpret it as 
    observation metadata, otherwise None.

    md_name is the column name for the last column if non-numeric

    NOTE: this is intended to be close to how QIIME classic OTU tables are
    parsed with the exception of the additional md_name field

    This function is ported from QIIME (http://www.qiime.org), previously named
    parse_classic_otu_table. QIIME is a GPL project, but we obtained permission
    from the authors of this function to port it to the BIOM Format project
    (and keep it under BIOM's BSD license).
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

class MetadataMap(dict):

    @classmethod
    def fromFile(cls, lines, strip_quotes=True, suppress_stripping=False,
                 header=None, process_fns=None):
        """Parse mapping file that relates samples or observations to metadata.

        Format: header line with fields
                optionally other comment lines starting with #
                tab-delimited fields

        process_fns: a dictionary of functions to apply to metadata categories.
         the keys should be the column headings, and the values should be
         functions which take a single value. For example, if the values in a
         column called "taxonomy" should be split on semi-colons before being
         added as metadata, and all other columns should be left as-is,
         process_fns should be:
          {'taxonomy': lambda x: x.split(';')}

        Assumes the first column in the mapping file is the id.

        This method is ported from QIIME (http://www.qiime.org), previously
        named parse_mapping_file/parse_mapping_file_to_dict. QIIME is a GPL
        project, but we obtained permission from the authors of this method
        to port it to the BIOM Format project (and keep it under BIOM's BSD
        license).
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
        
        # if the user didn't provide process functions, initialize as 
        # an empty dict
        if process_fns == None:
            process_fns = {}

        # Create lists to store the results
        mapping_data = []
        header = header or []
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
            raise BiomParseException("No header line was found in mapping "
                                     "file.")
        if not mapping_data:
            raise BiomParseException("No data found in mapping file.")

        first_col = [i[0] for i in mapping_data]
        if len(first_col) != len(set(first_col)):
            raise BiomParseException("First column values are not unique! "
                                     "Cannot be ids.")

        mapping = {}
        for vals in mapping_data:
            current_d = {}
            for k,v in zip(header[1:], vals[1:]):
                try:
                    current_d[k] = process_fns[k](v)
                except KeyError:
                    current_d[k] = v
            mapping[vals[0]] = current_d

        return cls(mapping)

    def __init__(self, mapping):
        """Accepts dictionary mapping IDs to metadata.

        ``mapping`` should be a dictionary mapping an ID to a dictionary of
        metadata. For example:

        {'Sample1': {'Treatment': 'Fast'}, 'Sample2': {'Treatment': 'Control'}}
        """
        super(MetadataMap, self).__init__(mapping)

def generatedby():
    """Returns a generated by string"""
    return 'BIOM-Format %s' % __version__

def parse_biom_table_hdf5(fp):
    table_f = h5py.File(fp, 'r')

    obs_ids = table_f['observations/ids'].value
    if 'observations/metadata' in table_f:
        obs_md = json.loads(table_f['observations/metadata'].value)

    sample_ids = table_f['samples/ids'].value
    if 'samples/metadata' in table_f:
        sample_md = json.loads(table_f['samples/metadata'].value)

    table_id = table_f.attrs['id']
    table_type = table_f.attrs['type']

    sparse_obj = SparseObj(*table_f.attrs['shape'])
    data_grp = table_f['data']
    sparse_obj._matrix = coo_matrix((data_grp['values'].value,
                                    (data_grp['rows'].value,
                                     data_grp['columns'].value)),
                                    shape=table_f.attrs['shape'])
    table_f.close()

    return Table(sparse_obj, sample_ids, obs_ids, sample_md, obs_md, 
                table_id)

def convert_table_to_biom(table_f,sample_mapping, obs_mapping, 
        process_func, **kwargs):
    """Convert a contigency table to a biom table
    
    sample_mapping : dict of {'sample_id':metadata} or None
    obs_mapping : dict of {'obs_id':metadata} or None
    process_func: a function to transform observation metadata
    dtype : type of table data
    """
    otu_table = parse_classic_table_to_rich_table(table_f, sample_mapping, 
                                                  obs_mapping, process_func,
                                                  **kwargs)
    return otu_table.getBiomFormatJsonString(generatedby())

def biom_meta_to_string(metadata, replace_str=':'):
    """ Determine which format the metadata is (e.g. str, list, or list of lists) and then convert to a string"""

    #Note that since ';' and '|' are used as seperators we must replace them if they exist
  
    #metadata is just a string (not a list)
    if isinstance(metadata,str) or isinstance(metadata,unicode):
        return metadata.replace(';',replace_str)

    elif isinstance(metadata,list):
        
        #metadata is list of lists
        if isinstance(metadata[0], list):
            new_metadata=[]
            for x in metadata:
                #replace erroneus delimiters
                values=[y.replace(';',replace_str).replace('|',replace_str).strip() for y in x]
                new_metadata.append("; ".join(values))
            return "|".join(new_metadata)

        #metadata is list (of strings)
        else:
            return "; ".join(x.replace(';',replace_str).strip() for x in metadata)

def convert_biom_to_table(biom_f, header_key=None, header_value=None, \
        md_format=None):
    """Convert a biom table to a contigency table"""
    table = parse_biom_table(biom_f)

    if md_format is None:
        md_format = biom_meta_to_string

    if table.ObservationMetadata is None:
        return table.delimitedSelf()
    
    if header_key in table.ObservationMetadata[0]:
        return table.delimitedSelf(header_key=header_key, 
                                       header_value=header_value,
                                       metadata_formatter=md_format)
    else:
        return table.delimitedSelf()
