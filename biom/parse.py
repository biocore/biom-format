#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division
import os
from string import maketrans
import numpy as np
from biom import __version__
from biom.exception import BiomParseException
from biom.table import table_factory, nparray_to_sparseobj
from functools import partial
import json
from numpy import asarray
from scipy.sparse import csr_matrix, csc_matrix
from biom.backends.scipysparse import ScipySparseMat

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Justin Kuczynski", "Daniel McDonald", "Greg Caporaso",
               "Jose Carlos Clemente Litran", "Adam Robbins-Pianka"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

MATRIX_ELEMENT_TYPE = {'int': int, 'float': float, 'unicode': unicode,
                       u'int': int, u'float': float, u'unicode': unicode}

QUOTE = '"'
JSON_OPEN = set(["[", "{"])
JSON_CLOSE = set(["]", "}"])
JSON_SKIP = set([" ", "\t", "\n", ","])
JSON_START = set(
    ["0",
     "1",
     "2",
     "3",
     "4",
     "5",
     "6",
     "7",
     "8",
     "9",
     "{",
     "[",
     '"'])


def direct_parse_key(biom_str, key):
    """Returns key:value from the biom string, or ""

    This method pulls an arbitrary key/value pair out from a BIOM string
    """
    base_idx = biom_str.find('"%s":' % key)
    if base_idx == -1:
        return ""
    else:
        start_idx = base_idx + len(key) + 3  # shift over "key":

    # find the start token
    cur_idx = start_idx
    while biom_str[cur_idx] not in JSON_START:
        cur_idx += 1

    if biom_str[cur_idx] not in JSON_OPEN:
        # do we have a number?
        while biom_str[cur_idx] not in [",", "{", "}"]:
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
                except IndexError:  # got an int or float?
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
    if axis not in ['observations', 'samples']:
        raise IndexError("Unknown axis type")

    # it would be nice if all of these lookups could be done in a single
    # traversal of biom_str, but it likely is at the cost of code complexity
    shape_kv_pair = direct_parse_key(biom_str, "shape")
    if shape_kv_pair == "":
        raise ValueError("biom_str does not appear to be in BIOM format!")

    data_fields = direct_parse_key(biom_str, "data")
    if data_fields == "":
        raise ValueError("biom_str does not appear to be in BIOM format!")

    matrix_type_kv_pair = direct_parse_key(biom_str, "matrix_type")
    if matrix_type_kv_pair == "":
        raise ValueError("biom_str does not appear to be in BIOM format!")

    # determine shape
    raw_shape = shape_kv_pair.split(':')[-1].replace("[", "").replace("]", "")
    n_rows, n_cols = map(int, raw_shape.split(","))

    # slice to just data
    data_start = data_fields.find('[') + 1
    # trim trailing ]
    data_fields = data_fields[data_start:len(data_fields) - 1]

    # bounds check
    if min(to_keep) < 0:
        raise IndexError("Observations to keep are out of bounds!")

    # more bounds check and set new shape
    new_shape = "[%d, %d]"
    if axis == 'observations':
        if max(to_keep) >= n_rows:
            raise IndexError("Observations to keep are out of bounds!")
        new_shape = new_shape % (len(to_keep), n_cols)
    elif axis == 'samples':
        if max(to_keep) >= n_cols:
            raise IndexError("Samples to keep are out of bounds!")
        new_shape = new_shape % (n_rows, len(to_keep))

    to_keep = set(to_keep)
    new_data = []

    if axis == 'observations':
        new_data = _direct_slice_data_sparse_obs(data_fields, to_keep)
    elif axis == 'samples':
        new_data = _direct_slice_data_sparse_samp(data_fields, to_keep)

    return '"data": %s, "shape": %s' % (new_data, new_shape)

STRIP_F = lambda x: x.strip("[] \n\t")


def _remap_axis_sparse_obs(rcv, lookup):
    """Remap a sparse observation axis"""
    row, col, value = map(STRIP_F, rcv.split(','))
    return "%s,%s,%s" % (lookup[row], col, value)


def _remap_axis_sparse_samp(rcv, lookup):
    """Remap a sparse sample axis"""
    row, col, value = map(STRIP_F, rcv.split(','))
    return "%s,%s,%s" % (row, lookup[col], value)


def _direct_slice_data_sparse_obs(data, to_keep):
    """slice observations from data

    data : raw data string from a biom file
    to_keep : rows to keep
    """
    # interogate all the datas
    new_data = []
    remap_lookup = dict([(str(v), i) for i, v in enumerate(sorted(to_keep))])
    for rcv in data.split('],'):
        r, c, v = STRIP_F(rcv).split(',')
        if r in remap_lookup:
            new_data.append(_remap_axis_sparse_obs(rcv, remap_lookup))
    return '[[%s]]' % '],['.join(new_data)


def _direct_slice_data_sparse_samp(data, to_keep):
    """slice samples from data

    data : raw data string from a biom file
    to_keep : columns to keep
    """
    # could do sparse obs/samp in one forloop, but then theres the
    # expense of the additional if-statement in the loop
    new_data = []
    remap_lookup = dict([(str(v), i) for i, v in enumerate(sorted(to_keep))])
    for rcv in data.split('],'):
        r, c, v = rcv.split(',')
        if c in remap_lookup:
            new_data.append(_remap_axis_sparse_samp(rcv, remap_lookup))
    return '[[%s]]' % '],['.join(new_data)


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
        raise ValueError("Unknown axis!")

    if axis_data == "":
        raise ValueError("biom_str does not appear to be in BIOM format!")

    axis_data = json.loads("{%s}" % axis_data)

    all_ids = set([v['id'] for v in axis_data[axis_key]])
    if not to_keep.issubset(all_ids):
        raise KeyError("Not all of the to_keep ids are in biom_str!")

    idxs = [i for i, v in enumerate(axis_data[axis_key]) if v['id'] in to_keep]
    idxs_lookup = set(idxs)

    subset = {axis_key: []}
    for i, v in enumerate(axis_data[axis_key]):
        if i in idxs_lookup:
            subset[axis_key].append(v)

    return idxs, json.dumps(subset)[1:-1]  # trim off { and }


def parse_biom_table(fp):
    if hasattr(fp, 'read'):
        return parse_biom_table_json(json.load(fp))
    elif isinstance(fp, list):
        return parse_biom_table_json(json.loads(''.join(fp)))
    else:
        return parse_biom_table_json(json.loads(fp))


def parse_biom_table_hdf5(h5grp, order='observation'):
    """Parse an HDF5 formatted BIOM table

    The expected structure of this group is below. A few basic definitions,
    N is the number of observations and M is the number of samples. Data are
    stored in both compressed sparse row (for observation oriented operations)
    and compressed sparse column (for sample oriented operations).

    ### ADD IN SCIPY SPARSE CSC/CSR URLS
    ### ADD IN WIKIPEDIA PAGE LINK TO CSR
    ### ALL THESE INTS CAN BE UINT, SCIPY DOES NOT BY DEFAULT STORE AS THIS
    ###     THOUGH
    ### METADATA ARE NOT REPRESENTED HERE YET
    ./id                     : str, an arbitrary ID
    ./type                   : str, the table type (e.g, OTU table)
    ./format-url             : str, a URL that describes the format
    ./format-version         : two element tuple of int32, major and minor
    ./generated-by           : str, what generated this file
    ./creation-date          : str, ISO format
    ./shape                  : two element tuple of int32, N by M
    ./nnz                    : int32 or int64, number of non zero elements
    ./observation            : Group
    ./observation/ids        : (N,) dataset of str or vlen str
    ./observation/data       : (N,) dataset of float64
    ./observation/indices    : (N,) dataset of int32
    ./observation/indptr     : (M+1,) dataset of int32
    [./observation/metadata] : Optional, JSON str, in index order with ids
    ./sample                 : Group
    ./sample/ids             : (M,) dataset of str or vlen str
    ./sample/data            : (M,) dataset of float64
    ./sample/indices         : (M,) dataset of int32
    ./sample/indptr          : (N+1,) dataset of int32
    [./sample/metadata]      : Optional, JSON str, in index order with ids
    Paramters
    ---------
    h5grp : a h5py ``Group`` or an open h5py ``File``
    order : 'observation' or 'sample' to indicate which data ordering to load
        the table as

    Returns
    -------
    Table
        A BIOM ``Table`` object

    See Also
    --------
    Table.format_hdf5

    Examples
    --------
    ### is it okay to actually create files in doctest?

    """
    if order not in ('observation', 'sample'):
        raise ValueError("Unknown order %s!" % order)

    # fetch all of the IDs
    obs_ids = h5grp['observation/ids'][:]
    samp_ids = h5grp['sample/ids'][:]

    # fetch all of the metadata
    no_md = np.array(["[]"])
    obs_md = json.loads(h5grp['observation'].get('metadata', no_md)[0])
    samp_md = json.loads(h5grp['sample'].get('metadata', no_md)[0])

    # construct the sparse representation
    rep = ScipySparseMat(len(obs_ids), len(samp_ids))

    # load the data
    data_path = partial(os.path.join, order)
    data = h5grp[data_path("data")]
    indices = h5grp[data_path("indices")]
    indptr = h5grp[data_path("indptr")]
    cs = (data, indices, indptr)
    rep._matrix = csc_matrix(cs) if order == 'sample' else csr_matrix(cs)

    return table_factory(rep, samp_ids, obs_ids, samp_md or None,
                         obs_md or None)


def parse_biom_table_json(json_table, data_pump=None):
    """Parse a biom otu table type"""
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
    complex_metadata = []
    for y in x.split('|'):
        simple_metadata = []
        for e in y.split(';'):
            simple_metadata.append(e.strip())
        complex_metadata.append(simple_metadata)
    return complex_metadata


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
        obs_metadata = [{t_md_name: process_func(v)} for v in t_md]

    if sample_mapping is None:
        sample_metadata = None
    else:
        sample_metadata = [sample_mapping[sample_id]
                           for sample_id in sample_ids]

    # will override any metadata from parsed table
    if obs_mapping is not None:
        obs_metadata = [obs_mapping[obs_id] for obs_id in obs_ids]

    data = nparray_to_sparseobj(data)

    return table_factory(data, sample_ids, obs_ids, sample_metadata,
                         obs_metadata)


def parse_classic_table(lines, delim='\t', dtype=float, header_mark=None,
                        md_parse=None):
    """Parse a classic table into (sample_ids, obs_ids, data, metadata, name)

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
            raise BiomParseException(
                "Input needs to support readlines or be indexable")

    # find header, the first line that is not empty and does not start with a #
    for idx, l in enumerate(lines):
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
            header = lines[idx - 1].strip().split(delim)[1:]

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
    def from_file(cls, lines, strip_quotes=True, suppress_stripping=False,
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
        if hasattr(lines, "upper"):
            # Try opening if a string was passed
            try:
                lines = open(lines, 'U')
            except IOError:
                raise BiomParseException("A string was passed that doesn't "
                                         "refer to an accessible filepath.")

        if strip_quotes:
            if suppress_stripping:
                # remove quotes but not spaces
                strip_f = lambda x: x.replace('"', '')
            else:
                # remove quotes and spaces
                strip_f = lambda x: x.replace('"', '').strip()
        else:
            if suppress_stripping:
                # don't remove quotes or spaces
                strip_f = lambda x: x
            else:
                # remove spaces but not quotes
                strip_f = lambda x: x.strip()

        # if the user didn't provide process functions, initialize as
        # an empty dict
        if process_fns is None:
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
                if len(tmp_line) < len(header):
                    tmp_line.extend([''] * (len(header) - len(tmp_line)))
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
            for k, v in zip(header[1:], vals[1:]):
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


def convert_table_to_biom(table_f, sample_mapping, obs_mapping,
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
    return otu_table.get_biom_format_json_string(generatedby())


def biom_meta_to_string(metadata, replace_str=':'):
    """Determine which format the metadata is and then convert to a string"""

    # Note that since ';' and '|' are used as seperators we must replace them
    # if they exist

    # metadata is just a string (not a list)
    if isinstance(metadata, str) or isinstance(metadata, unicode):
        return metadata.replace(';', replace_str)

    elif isinstance(metadata, list):
        transtab = maketrans(';|', replace_str)
        # metadata is list of lists
        if isinstance(metadata[0], list):
            new_metadata = []
            for x in metadata:
                # replace erroneus delimiters
                values = [y.strip().trans(transtab) for y in x]
                new_metadata.append("; ".join(values))
            return "|".join(new_metadata)

        # metadata is list (of strings)
        else:
            return (
                "; ".join(x.replace(';', replace_str).strip()
                          for x in metadata)
            )


def convert_biom_to_table(biom_f, header_key=None, header_value=None,
                          md_format=None):
    """Convert a biom table to a contigency table"""
    table = parse_biom_table(biom_f)

    if md_format is None:
        md_format = biom_meta_to_string

    if table.observation_metadata is None:
        return table.delimited_self()

    if header_key in table.observation_metadata[0]:
        return table.delimited_self(header_key=header_key,
                                    header_value=header_value,
                                    metadata_formatter=md_format)
    else:
        return table.delimited_self()
