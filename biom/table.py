#!/usr/bin/env python
"""
BIOM Table (:mod:`biom.table`)
==============================

The biom-format project provides rich ``Table`` objects to support use of the
BIOM file format. The objects encapsulate matrix data (such as OTU counts) and
abstract the interaction away from the programmer.

.. currentmodule:: biom.table

Classes
-------

.. autosummary::
   :toctree: generated/

   Table

Examples
--------
First, let's create a toy table to play around with. For this example, we're
going to construct a 10x4 `Table`, or one that has 10 observations and 4
samples. Each observation and sample will be given an arbitrary but unique
name. We'll also add on some metadata.

>>> import numpy as np
>>> from biom.table import Table
>>> data = np.arange(40).reshape(10, 4)
>>> sample_ids = ['S%d' % i for i in range(4)]
>>> observ_ids = ['O%d' % i for i in range(10)]
>>> sample_metadata = [{'environment': 'A'}, {'environment': 'B'},
...                    {'environment': 'A'}, {'environment': 'B'}]
>>> observ_metadata = [{'taxonomy': ['Bacteria', 'Firmicutes']},
...                    {'taxonomy': ['Bacteria', 'Firmicutes']},
...                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
...                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
...                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
...                    {'taxonomy': ['Bacteria', 'Bacteroidetes']},
...                    {'taxonomy': ['Bacteria', 'Bacteroidetes']},
...                    {'taxonomy': ['Bacteria', 'Firmicutes']},
...                    {'taxonomy': ['Bacteria', 'Firmicutes']},
...                    {'taxonomy': ['Bacteria', 'Firmicutes']}]
>>> table = Table(data, observ_ids, sample_ids, observ_metadata,
...               sample_metadata, table_id='Example Table')

Now that we have a table, let's explore it at a high level first.

>>> table
10 x 4 <class 'biom.table.Table'> with 39 nonzero entries (97% dense)
>>> print(table) # doctest: +NORMALIZE_WHITESPACE
# Constructed from biom file
#OTU ID S0  S1  S2  S3
O0  0.0 1.0 2.0 3.0
O1  4.0 5.0 6.0 7.0
O2  8.0 9.0 10.0    11.0
O3  12.0    13.0    14.0    15.0
O4  16.0    17.0    18.0    19.0
O5  20.0    21.0    22.0    23.0
O6  24.0    25.0    26.0    27.0
O7  28.0    29.0    30.0    31.0
O8  32.0    33.0    34.0    35.0
O9  36.0    37.0    38.0    39.0
>>> print(table.ids()) # doctest: +NORMALIZE_WHITESPACE
['S0' 'S1' 'S2' 'S3']
>>> print(table.ids(axis='observation')) # doctest: +NORMALIZE_WHITESPACE
['O0' 'O1' 'O2' 'O3' 'O4' 'O5' 'O6' 'O7' 'O8' 'O9']
>>> print(table.nnz)  # number of nonzero entries
39

While it's fun to just poke at the table, let's dig deeper. First, we're going
to convert `table` into relative abundances (within each sample), and then
filter `table` to just the samples associated with environment 'A'. The
filtering gets fancy: we can pass in an arbitrary function to determine what
samples we want to keep. This function must accept a sparse vector of values,
the corresponding ID and the corresponding metadata, and should return ``True``
or ``False``, where ``True`` indicates that the vector should be retained.

>>> normed = table.norm(axis='sample', inplace=False)
>>> filter_f = lambda values, id_, md: md['environment'] == 'A'
>>> env_a = normed.filter(filter_f, axis='sample', inplace=False)
>>> print(env_a) # doctest: +NORMALIZE_WHITESPACE
# Constructed from biom file
#OTU ID S0  S2
O0  0.0 0.01
O1  0.0222222222222 0.03
O2  0.0444444444444 0.05
O3  0.0666666666667 0.07
O4  0.0888888888889 0.09
O5  0.111111111111  0.11
O6  0.133333333333  0.13
O7  0.155555555556  0.15
O8  0.177777777778  0.17
O9  0.2 0.19

But, what if we wanted individual tables per environment? While we could just
perform some fancy iteration, we can instead just rely on `Table.partition` for
these operations. `partition`, like `filter`, accepts a function. However, the
`partition` method only passes the corresponding ID and metadata to the
function. The function should return what partition the data are a part of.
Within this example, we're also going to sum up our tables over the partitioned
samples. Please note that we're using the original table (ie, not normalized)
here.

>>> part_f = lambda id_, md: md['environment']
>>> env_tables = table.partition(part_f, axis='sample')
>>> for partition, env_table in env_tables:
...     print(partition, env_table.sum('sample'))
A [ 180.  200.]
B [ 190.  210.]

For this last example, and to highlight a bit more functionality, we're going
to first transform the table such that all multiples of three will be retained,
while all non-multiples of three will get set to zero. Following this, we'll
then collpase the table by taxonomy, and then convert the table into
presence/absence data.

First, let's setup the transform. We're going to define a function that takes
the modulus of every value in the vector, and see if it is equal to zero. If it
is equal to zero, we'll keep the value, otherwise we'll set the value to zero.

>>> transform_f = lambda v,i,m: np.where(v % 3 == 0, v, 0)
>>> mult_of_three = tform = table.transform(transform_f, inplace=False)
>>> print(mult_of_three) # doctest: +NORMALIZE_WHITESPACE
# Constructed from biom file
#OTU ID S0  S1  S2  S3
O0  0.0 0.0 0.0 3.0
O1  0.0 0.0 6.0 0.0
O2  0.0 9.0 0.0 0.0
O3  12.0    0.0 0.0 15.0
O4  0.0 0.0 18.0    0.0
O5  0.0 21.0    0.0 0.0
O6  24.0    0.0 0.0 27.0
O7  0.0 0.0 30.0    0.0
O8  0.0 33.0    0.0 0.0
O9  36.0    0.0 0.0 39.0

Next, we're going to collapse the table over the phylum level taxon. To do
this, we're going to define a helper variable for the index position of the
phylum (see the construction of the table above). Next, we're going to pass
this to `Table.collapse`, and since we want to collapse over the observations,
we'll need to specify 'observation' as the axis.

>>> phylum_idx = 1
>>> collapse_f = lambda id_, md: '; '.join(md['taxonomy'][:phylum_idx + 1])
>>> collapsed = mult_of_three.collapse(collapse_f, axis='observation')
>>> print(collapsed) # doctest: +NORMALIZE_WHITESPACE
# Constructed from biom file
#OTU ID S0  S1  S2  S3
Bacteria; Firmicutes  7.2 6.6 7.2 8.4
Bacteria; Proteobacteria  4.0 3.0 6.0 5.0
Bacteria; Bacteroidetes   12.0    10.5    0.0 13.5

Finally, let's convert the table to presence/absence data.

>>> pa = collapsed.pa()
>>> print(pa) # doctest: +NORMALIZE_WHITESPACE
# Constructed from biom file
#OTU ID S0  S1  S2  S3
Bacteria; Firmicutes  1.0 1.0 1.0 1.0
Bacteria; Proteobacteria  1.0 1.0 1.0 1.0
Bacteria; Bacteroidetes   1.0 1.0 0.0 1.0

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2020, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np
import scipy.stats
import h5py
from copy import deepcopy
from datetime import datetime
from json import dumps as _json_dumps, JSONEncoder
from functools import reduce, partial
from operator import itemgetter
from collections import defaultdict
from collections.abc import Hashable, Iterable
from numpy import ndarray, asarray, zeros, newaxis
from scipy.sparse import (coo_matrix, csc_matrix, csr_matrix, isspmatrix,
                          vstack, hstack, dok_matrix)
import pandas as pd
import re
from biom.exception import (TableException, UnknownAxisError, UnknownIDError,
                            DisjointIDError)
from biom.util import (get_biom_format_version_string,
                       get_biom_format_url_string, flatten, natsort,
                       prefer_self, index_list, H5PY_VLEN_STR,
                       __format_version__)
from biom.err import errcheck
from ._filter import _filter
from ._transform import _transform
from ._subsample import subsample


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2020, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Adam Robbins-Pianka",
               "Joshua Shorenstein", "Jose Antonio Navas Molina",
               "Jorge Cañardo Alastuey", "Steven Brown"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


MATRIX_ELEMENT_TYPE = {'int': int, 'float': float, 'unicode': str,
                       'int': int, 'float': float, 'unicode': str}


# NpEncoder from:
# https://stackoverflow.com/a/57915246/19741
class NpEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


dumps = partial(_json_dumps, cls=NpEncoder)


def _identify_bad_value(dtype, fields):
    """Identify the first value which cannot be cast

    Paramters
    ---------
    dtype : type
        A type to cast to
    fields : Iterable of str
        A series of str to cast into dtype

    Returns
    -------
    str or None
        A value that cannot be cast
    int or None
        The index of the value that cannot be cast
    """
    badval = None
    badidx = None
    for idx, v in enumerate(fields):
        try:
            dtype(v)
        except:  # noqa
            badval = v
            badidx = idx
            break
    return (badval, badidx)


def general_parser(x):
    if isinstance(x, bytes):
        x = x.decode('utf8')
    return x


def vlen_list_of_str_parser(value):
    """Parses the taxonomy value"""
    new_value = []
    for v in value:
        if v:
            if isinstance(v, bytes):
                v = v.decode('utf8')
            new_value.append(v)

    return new_value if new_value else None


def general_formatter(grp, header, md, compression):
    """Creates a dataset for a general atomic type category"""
    shape = (len(md),)
    dtypes = [type(m[header]) for m in md]

    # "/" are considered part of the path in hdf5 and must be
    # escaped. However, escaping with "\" leads to a truncation
    # so let's replace "/" with an unexpected keyword.
    sanitized = header.replace('/', '@@SLASH@@')
    name = 'metadata/%s' % sanitized

    if set(dtypes).issubset({str}):
        grp.create_dataset(name, shape=shape,
                           dtype=H5PY_VLEN_STR,
                           data=[m[header].encode('utf8') for m in md],
                           compression=compression)
    elif set(dtypes).issubset({list, tuple}):
        vlen_list_of_str_formatter(grp, header, md, compression)
    else:
        formatted = []
        dtypes_used = []
        for dt, m in zip(dtypes, md):
            val = m[header]
            if val is None:
                val = ''
                dt = str

            if dt == str:
                val = val.encode('utf8')

            formatted.append(val)
            dtypes_used.append(dt)

        if set(dtypes_used).issubset({str}):
            dtype_to_use = H5PY_VLEN_STR
        else:
            dtype_to_use = None

        # try our best...
        grp.create_dataset(
            name, shape=(len(md),),
            dtype=dtype_to_use,
            data=formatted,
            compression=compression)


def vlen_list_of_str_formatter(grp, header, md, compression):
    """Creates a (N, ?) vlen str dataset"""
    # It is possible that the value for some sample/observation
    # is None. In that case, we still need to see them as
    # iterables, but their length will be 0
    iterable_checks = []
    lengths = []
    for m in md:
        if m[header] is None:
            iterable_checks.append(True)
        elif isinstance(m.get(header), str):
            iterable_checks.append(False)
        else:
            iterable_checks.append(
                isinstance(m.get(header, []), Iterable))
            lengths.append(len(m[header]))

    if not np.all(iterable_checks):
        if header == 'taxonomy':
            # attempt to handle the general case issue where the taxonomy
            # was not split on semicolons and represented as a flat string
            # instead of a list
            def split_and_strip(i):
                parts = i.split(';')
                return [p.strip() for p in parts]
            try:
                new_md = []
                lengths = []
                for m in md:
                    parts = split_and_strip(m[header])
                    new_md.append({header: parts})
                    lengths.append(len(parts))
                md = new_md
            except:  # noqa
                raise TypeError("Category '%s' is not formatted properly. The "
                                "most common issue is when 'taxonomy' is "
                                "represented as a flat string instead of a "
                                "list. An attempt was made to split this "
                                "field on a ';' to coerce it into a list but "
                                "it failed. An example entry (which is not "
                                "assured to be the problematic entry) is "
                                "below:\n%s" % (header, md[0][header]))
        else:
            raise TypeError(
                "Category %s not formatted correctly. Did you pass"
                " --process-obs-metadata taxonomy when converting "
                " from tsv? Please see Table.to_hdf5 docstring for"
                " more information" % header)

    max_list_len = max(lengths)
    shape = (len(md), max_list_len)
    data = np.empty(shape, dtype=object)
    for i, m in enumerate(md):
        if m[header] is None:
            continue
        value = np.asarray(m[header])
        data[i, :len(value)] = [v.encode('utf8') for v in value]
    # Change the None entries on data to empty strings ""
    data = np.where(data == np.array(None), "", data)
    grp.create_dataset(
        'metadata/%s' % header, shape=shape,
        dtype=H5PY_VLEN_STR, data=data,
        compression=compression)


class Table:
    """The (canonically pronounced 'teh') Table.

    Give in to the power of the Table!

    Creates an in-memory representation of a BIOM file. BIOM version 1.0 is
    based on JSON to provide the overall structure for the format while
    versions 2.0 and 2.1 are based on HDF5. For more information see [1]_
    and [2]_

    Paramaters
    ----------
    data : array_like
        An (N,M) sample by observation matrix represented as one of these
        types:
        * An 1-dimensional array of values
        * An n-dimensional array of values
        * An empty list
        * A list of numpy arrays
        * A list of dict
        * A list of sparse matrices
        * A dictionary of values
        * A list of lists
        * A sparse matrix of values
    observation_ids : array_like of str
        A (N,) dataset of the observation IDs, where N is the total number
        of IDs
    sample_ids : array_like of str
        A (M,) dataset of the sample IDs, where M is the total number of IDs
    observation_metadata : list of dicts, optional
        per observation dictionary of annotations where every key represents a
        metadata field that contains specific metadata information,
        ie taxonomy, KEGG pathway, etc
    sample_metadata : array_like of dicts, optional
        per sample dictionary of annotations where every key represents a
        metadata field that contains sample specific metadata information, ie
    table_id : str, optional
        A field that can be used to identify the table
    type : str, see notes
        The type of table represented
    create_date : str, optional
        Date that this table was built
    generated_by : str, optional
        Individual who built the table
    observation_group_metadata : list, optional
        group that contains observation specific group metadata information
        (e.g., phylogenetic tree)
    sample_group_metadata : list, optional
        group that contains sample specific group metadata information
        (e.g., relationships between samples)

    Attributes
    ----------
    shape
    dtype
    nnz
    matrix_data
    type
    table_id
    create_date
    generated_by
    format_version

    Notes
    -----
    Allowed table types are None, "OTU table", "Pathway table", "Function
    table", "Ortholog table", "Gene table", "Metabolite table", "Taxon table"

    Raises
    ------
    TableException
        When an invalid table type is provided.

    References
    ----------
    .. [1] http://biom-format.org/documentation/biom_format.html
    .. [2] D. McDonald, et al. "The Biological Observation Matrix (BIOM) format
       or: how I learned to stop worrying and love the ome-ome"
       GigaScience 2012 1:7
    """

    def __init__(self, data, observation_ids, sample_ids,
                 observation_metadata=None, sample_metadata=None,
                 table_id=None, type=None, create_date=None, generated_by=None,
                 observation_group_metadata=None, sample_group_metadata=None,
                 validate=True, observation_index=None, sample_index=None,
                 **kwargs):

        self.type = type
        self.table_id = table_id
        self.create_date = create_date
        self.generated_by = generated_by
        self.format_version = __format_version__

        if not isspmatrix(data):
            shape = (len(observation_ids), len(sample_ids))
            input_is_dense = kwargs.get('input_is_dense', False)
            self._data = Table._to_sparse(data, input_is_dense=input_is_dense,
                                          shape=shape)
        else:
            self._data = data.tocsr()

        self._data = self._data.astype(float)

        self._sample_ids = np.asarray(sample_ids)
        self._observation_ids = np.asarray(observation_ids)

        if sample_metadata is not None:
            # not m will evaluate True if the object tested is None or
            # an empty dict, etc.
            if {not m for m in sample_metadata} == {True, }:
                self._sample_metadata = None
            else:
                self._sample_metadata = tuple(sample_metadata)
        else:
            self._sample_metadata = None

        if observation_metadata is not None:
            # not m will evaluate True if the object tested is None or
            # an empty dict, etc.
            if {not m for m in observation_metadata} == {True, }:
                self._observation_metadata = None
            else:
                self._observation_metadata = tuple(observation_metadata)
        else:
            self._observation_metadata = None

        self._sample_group_metadata = sample_group_metadata
        self._observation_group_metadata = observation_group_metadata

        if validate:
            errcheck(self)

        # These will be set by _index_ids()
        self._sample_index = None
        self._obs_index = None

        self._cast_metadata()
        self._index_ids(observation_index, sample_index)

    def _index_ids(self, observation_index, sample_index):
        """Sets lookups {id:index in _data}.

        Should only be called in constructor as this modifies state.
        """
        if sample_index is None:
            self._sample_index = index_list(self._sample_ids)
        else:
            self._sample_index = sample_index

        if observation_index is None:
            self._obs_index = index_list(self._observation_ids)
        else:
            self._obs_index = observation_index

    def _index(self, axis='sample'):
        """Return the index lookups of the given axis

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            Axis to get the index dict. Defaults to 'sample'

        Returns
        -------
        dict
            lookups {id:index}

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        """
        if axis == 'sample':
            return self._sample_index
        elif axis == 'observation':
            return self._obs_index
        else:
            raise UnknownAxisError(axis)

    def _conv_to_self_type(self, vals, transpose=False, dtype=None):
        """For converting vectors to a compatible self type"""
        if dtype is None:
            dtype = self.dtype

        if isspmatrix(vals):
            return vals
        else:
            return Table._to_sparse(vals, transpose, dtype)

    @staticmethod
    def _to_dense(vec):
        """Converts a row/col vector to a dense numpy array.

        Always returns a 1-D row vector for consistency with numpy iteration
        over arrays.
        """
        dense_vec = vec.toarray()

        if vec.shape == (1, 1):
            # Handle the special case where we only have a single element, but
            # we don't want to return a numpy scalar / 0-d array. We still want
            # to return a vector of length 1.
            return dense_vec.reshape(1)
        else:
            return np.squeeze(dense_vec)

    @staticmethod
    def _to_sparse(values, transpose=False, dtype=float, input_is_dense=False,
                   shape=None):
        """Try to return a populated scipy.sparse matrix.

        NOTE: assumes the max value observed in row and col defines the size of
        the matrix.
        """
        # if it is a vector
        if isinstance(values, ndarray) and len(values.shape) == 1:
            if transpose:
                mat = nparray_to_sparse(values[:, newaxis], dtype)
            else:
                mat = nparray_to_sparse(values, dtype)
            return mat
        if isinstance(values, ndarray):
            if transpose:
                mat = nparray_to_sparse(values.T, dtype)
            else:
                mat = nparray_to_sparse(values, dtype)
            return mat
        # the empty list
        elif isinstance(values, list) and len(values) == 0:
            return coo_matrix((0, 0))
        # list of np vectors
        elif isinstance(values, list) and isinstance(values[0], ndarray):
            mat = list_nparray_to_sparse(values, dtype)
            if transpose:
                mat = mat.T
            return mat
        # list of dicts, each representing a row in row order
        elif isinstance(values, list) and isinstance(values[0], dict):
            mat = list_dict_to_sparse(values, dtype)
            if transpose:
                mat = mat.T
            return mat
        # list of scipy.sparse matrices, each representing a row in row order
        elif isinstance(values, list) and isspmatrix(values[0]):
            mat = list_sparse_to_sparse(values, dtype)
            if transpose:
                mat = mat.T
            return mat
        elif isinstance(values, dict):
            mat = dict_to_sparse(values, dtype, shape)
            if transpose:
                mat = mat.T
            return mat
        elif isinstance(values, list) and isinstance(values[0], list):
            if input_is_dense:
                d = coo_matrix(values)
                mat = coo_arrays_to_sparse((d.data, (d.row, d.col)),
                                           dtype=dtype, shape=shape)
            else:
                mat = list_list_to_sparse(values, dtype, shape=shape)
            return mat
        elif isspmatrix(values):
            mat = values
            if transpose:
                mat = mat.transpose()
            return mat
        else:
            raise TableException("Unknown input type")

    def _cast_metadata(self):
        """Casts all metadata to defaultdict to support default values.

        Should be called after any modifications to sample/observation
        metadata.
        """
        def cast_metadata(md):
            """Do the actual casting"""
            default_md = []
            # if we have a list of [None], set to None
            if md is not None:
                if md.count(None) == len(md):
                    return None
            if md is not None:
                for item in md:
                    d = defaultdict(lambda: None)

                    if isinstance(item, dict):
                        d.update(item)
                    elif item is None:
                        pass
                    else:
                        raise TableException("Unable to cast metadata: %s" %
                                             repr(item))
                    default_md.append(d)
                return tuple(default_md)
            return md

        self._sample_metadata = cast_metadata(self._sample_metadata)
        self._observation_metadata = cast_metadata(self._observation_metadata)

        self._sample_group_metadata = (
            self._sample_group_metadata
            if self._sample_group_metadata else None)
        self._observation_group_metadata = (
            self._observation_group_metadata
            if self._observation_group_metadata else None)

    @property
    def shape(self):
        """The shape of the underlying contingency matrix"""
        return self._data.shape

    @property
    def dtype(self):
        """The type of the objects in the underlying contingency matrix"""
        return self._data.dtype

    @property
    def nnz(self):
        """Number of non-zero elements of the underlying contingency matrix"""
        self._data.eliminate_zeros()
        return self._data.nnz

    @property
    def matrix_data(self):
        """The sparse matrix object"""
        return self._data

    def length(self, axis='sample'):
        """Return the length of an axis

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            The axis to operate on

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> from biom import example_table
        >>> print(example_table.length(axis='sample'))
        3
        >>> print(example_table.length(axis='observation'))
        2
        """
        if axis not in ('sample', 'observation'):
            raise UnknownAxisError(axis)

        return self.shape[1] if axis == 'sample' else self.shape[0]

    def add_group_metadata(self, group_md, axis='sample'):
        """Take a dict of group metadata and add it to an axis

        Parameters
        ----------
        group_md : dict of tuples
            `group_md` should be of the form ``{category: (data type, value)``
        axis : {'sample', 'observation'}, optional
            The axis to operate on

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        """
        if axis == 'sample':
            if self._sample_group_metadata is not None:
                self._sample_group_metadata.update(group_md)
            else:
                self._sample_group_metadata = group_md
        elif axis == 'observation':
            if self._observation_group_metadata is not None:
                self._observation_group_metadata.update(group_md)
            else:
                self._observation_group_metadata = group_md
        else:
            raise UnknownAxisError(axis)

    def del_metadata(self, keys=None, axis='whole'):
        """Remove metadata from an axis

        Parameters
        ----------
        keys : list of str, optional
            The keys to remove from metadata. If None, all keys from the axis
            are removed.
        axis : {'sample', 'observation', 'whole'}, optional
            The axis to operate on. If 'whole', the operation is applied to
            both the sample and observation axes.

        Raises
        ------
        UnknownAxisError
            If the requested axis does not exist.

        Examples
        --------
        >>> from biom import Table
        >>> import numpy as np
        >>> tab = Table(np.array([[1, 2], [3, 4]]),
        ...             ['O1', 'O2'],
        ...             ['S1', 'S2'],
        ...             sample_metadata=[{'barcode': 'ATGC', 'env': 'A'},
        ...                              {'barcode': 'GGTT', 'env': 'B'}])
        >>> tab.del_metadata(keys=['env'])
        >>> for id, md in zip(tab.ids(), tab.metadata()):
        ...     print(id, list(md.items()))
        S1 [('barcode', 'ATGC')]
        S2 [('barcode', 'GGTT')]
        """
        if axis == 'whole':
            axes = ['sample', 'observation']
        elif axis in ('sample', 'observation'):
            axes = [axis]
        else:
            raise UnknownAxisError("%s is not recognized" % axis)

        if keys is None:
            if axis == 'whole':
                self._sample_metadata = None
                self._observation_metadata = None
            elif axis == 'sample':
                self._sample_metadata = None
            else:
                self._observation_metadata = None
            return

        for ax in axes:
            if self.metadata(axis=ax) is None:
                continue

            for i, md in zip(self.ids(axis=ax), self.metadata(axis=ax)):
                for k in keys:
                    if k in md:
                        del md[k]

            # for consistency with init on absence of metadata
            empties = {True if not md else False
                       for md in self.metadata(axis=ax)}
            if empties == {True, }:
                if ax == 'sample':
                    self._sample_metadata = None
                else:
                    self._observation_metadata = None

    def add_metadata(self, md, axis='sample'):
        """Take a dict of metadata and add it to an axis.

        Parameters
        ----------
        md : dict of dict
            `md` should be of the form ``{id: {dict_of_metadata}}``
        axis : {'sample', 'observation'}, optional
            The axis to operate on
        """
        metadata = self.metadata(axis=axis)
        if metadata is not None:
            for id_, md_entry in md.items():
                if self.exists(id_, axis=axis):
                    idx = self.index(id_, axis=axis)
                    metadata[idx].update(md_entry)
        else:
            ids = self.ids(axis=axis)
            if axis == 'sample':
                self._sample_metadata = tuple(
                    md[id_] if id_ in md else None for id_ in ids)
            elif axis == 'observation':
                self._observation_metadata = tuple(
                    md[id_] if id_ in md else None for id_ in ids)
            else:
                raise UnknownAxisError(axis)
        self._cast_metadata()

    def __getitem__(self, args):
        """Handles row or column slices

        Slicing over an individual axis is supported, but slicing over both
        axes at the same time is not supported. Partial slices, such as
        `foo[0, 5:10]` are not supported, however full slices are supported,
        such as `foo[0, :]`.

        Parameters
        ----------
        args : tuple or slice
            The specific element (by index position) to return or an entire
            row or column of the data.

        Returns
        -------
        float or spmatrix
            A float is return if a specific element is specified, otherwise a
            spmatrix object representing a vector of sparse data is returned.

        Raises
        ------
        IndexError
            - If the matrix is empty
            - If the arguments do not appear to be a tuple
            - If a slice on row and column is specified
            - If a partial slice is specified

        Notes
        -----
        Switching between slicing rows and columns is inefficient.  Slicing of
        rows requires a CSR representation, while slicing of columns requires a
        CSC representation, and transforms are performed on the data if the
        data are not in the required representation. These transforms can be
        expensive if done frequently.

        .. shownumpydoc
        """
        if self.is_empty():
            raise IndexError("Cannot retrieve an element from an empty/null "
                             "table.")

        try:
            row, col = args
        except:  # noqa
            raise IndexError("Must specify (row, col).")

        if isinstance(row, slice) and isinstance(col, slice):
            raise IndexError("Can only slice a single axis.")

        if isinstance(row, slice):
            if row.start is None and row.stop is None:
                return self._get_col(col)
            else:
                raise IndexError("Can only handle full : slices per axis.")
        elif isinstance(col, slice):
            if col.start is None and col.stop is None:
                return self._get_row(row)
            else:
                raise IndexError("Can only handle full : slices per axis.")
        else:
            if self._data.getformat() == 'coo':
                self._data = self._data.tocsr()

            return self._data[row, col]

    def _get_row(self, row_idx):
        """Return the row at ``row_idx``.

        A row vector will be returned as a scipy.sparse matrix in csr format.

        Notes
        -----
        Switching between slicing rows and columns is inefficient.  Slicing of
        rows requires a CSR representation, while slicing of columns requires a
        CSC representation, and transforms are performed on the data if the
        data are not in the required representation. These transforms can be
        expensive if done frequently.

        """
        self._data = self._data.tocsr()
        return self._data.getrow(row_idx)

    def _get_col(self, col_idx):
        """Return the column at ``col_idx``.

        A column vector will be returned as a scipy.sparse matrix in csc
        format.

        Notes
        -----
        Switching between slicing rows and columns is inefficient.  Slicing of
        rows requires a CSR representation, while slicing of columns requires a
        CSC representation, and transforms are performed on the data if the
        data are not in the required representation. These transforms can be
        expensive if done frequently.

        """
        self._data = self._data.tocsc()
        return self._data.getcol(col_idx)

    def align_to_dataframe(self, metadata, axis='sample'):
        """ Aligns dataframe against biom table, only keeping common ids.

        Parameters
        ----------
        metadata : pd.DataFrame
            The metadata, either respect to the sample metadata
            or observation metadata.
        axis : {'sample', 'observation'}
            The axis on which to operate.

        Returns
        -------
        biom.Table
            A filtered biom table.
        pd.DataFrame
            A filtered metadata table.

        Examples
        --------
        >>> from biom import Table
        >>> import numpy as np
        >>> import pandas as pd
        >>> table = Table(np.array([[0, 0, 1, 1],
        ...                         [2, 2, 4, 4],
        ...                         [5, 5, 3, 3],
        ...                         [0, 0, 0, 1]]),
        ...               ['o1', 'o2', 'o3', 'o4'],
        ...               ['s1', 's2', 's3', 's4'])
        >>> metadata = pd.DataFrame([['a', 'control'],
        ...                          ['c', 'diseased'],
        ...                          ['b', 'control']],
        ...                         index=['s1', 's3', 's2'],
        ...                         columns=['Barcode', 'Treatment'])
        >>> res_table, res_metadata = table.align_to_dataframe(metadata)
        >>> print(res_table)
        # Constructed from biom file
        #OTU ID	s1	s2	s3
        o1	0.0	0.0	1.0
        o2	2.0	2.0	4.0
        o3	5.0	5.0	3.0
        >>> print(res_metadata)
           Barcode Treatment
        s1       a   control
        s2       b   control
        s3       c  diseased
        """
        ids = set(self.ids(axis=axis)) & set(metadata.index)
        if len(ids) == 0:
            raise TableException("No common ids between table and dataframe.")

        t = self.filter(ids, axis=axis, inplace=False)
        t.remove_empty()
        md = metadata.loc[t.ids(axis=axis)]
        return t, md

    def align_tree(self, tree, axis='observation'):
        r""" Aligns biom table against tree, only keeping common ids.

        Parameters
        ----------
        tree : skbio.TreeNode
            The tree object, either respect to the sample metadata
            or observation metadata.
        axis : {'sample', 'observation'}
            The axis on which to operate.

        Returns
        -------
        biom.Table
            A filtered biom table.
        skbio.TreeNode
            A filtered skbio TreeNode object.

        Examples
        --------
        >>> from biom import Table
        >>> import numpy as np
        >>> from skbio import TreeNode
        >>> table = Table(np.array([[0, 0, 1, 1],
        ...                         [2, 2, 4, 4],
        ...                         [5, 5, 3, 3],
        ...                         [0, 0, 0, 1]]),
        ...               ['o1', 'o2', 'o3', 'o4'],
        ...               ['s1', 's2', 's3', 's4'])
        >>> tree = TreeNode.read([u"((o1,o2)f,o3)r;"])
        >>> res_table, res_tree = table.align_tree(tree)
        >>> print(res_table)
        # Constructed from biom file
        #OTU ID	s1	s2	s3	s4
        o1	0.0	0.0	1.0	1.0
        o2	2.0	2.0	4.0	4.0
        o3	5.0	5.0	3.0	3.0
        >>> print(res_tree.ascii_art())
                            /-o1
                  /f-------|
        -r-------|          \-o2
                 |
                  \-o3
        """
        tips = {x.name for x in tree.tips()}
        common_tips = tips & set(self.ids(axis=axis))
        if len(common_tips) == 0:
            raise TableException("No common ids between table and tree.")
        _tree = tree.shear(names=common_tips)
        _table = self.filter(common_tips, axis=axis, inplace=False)
        _tree.prune()
        order = [n.name for n in _tree.tips()]
        _table = _table.sort_order(order, axis=axis)
        return _table, _tree

    def reduce(self, f, axis):
        """Reduce over axis using function `f`

        Parameters
        ----------
        f : function
            The function to use for the reduce operation
        axis : {'sample', 'observation'}
            The axis on which to operate

        Returns
        -------
        numpy.array
            A one-dimensional array representing the reduced rows
            (observations) or columns (samples) of the data matrix

        Raises
        ------
        UnknownAxisError
            If `axis` is neither "sample" nor "observation"
        TableException
            If the table's data matrix is empty

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 table

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               [{'foo': 'bar'}, {'x': 'y'}], None)

        Create a reduce function

        >>> func = lambda x, y: x + y

        Reduce table on samples

        >>> table.reduce(func, 'sample') # doctest: +NORMALIZE_WHITESPACE
        array([  1.,   3.,  43.])

        Reduce table on observations

        >>> table.reduce(func, 'observation') # doctest: +NORMALIZE_WHITESPACE
        array([  1.,  46.])
        """
        if self.is_empty():
            raise TableException("Cannot reduce an empty table")

        # np.apply_along_axis might reduce type conversions here and improve
        # speed. am opting for reduce right now as I think its more readable
        return asarray([reduce(f, v) for v in self.iter_data(axis=axis)])

    def sum(self, axis='whole'):
        """Returns the sum by axis

        Parameters
        ----------
        axis : {'whole', 'sample', 'observation'}, optional
            The axis on which to operate.

        Returns
        -------
        numpy.array or float
            If `axis` is "whole", returns an float representing the whole
            table sum. If `axis` is either "sample" or "observation", returns a
            numpy.array that holds a sum for each sample or observation,
            respectively.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Add all values in the table:

        >>> table.sum()
        array(47.0)

        Add all values per sample:

        >>> table.sum(axis='sample') # doctest: +NORMALIZE_WHITESPACE
        array([  1.,  3.,  43.])

        Add all values per observation:

        >>> table.sum(axis='observation') # doctest: +NORMALIZE_WHITESPACE
        array([  1.,  46.])
        """
        if axis == 'whole':
            axis = None
        elif axis == 'sample':
            axis = 0
        elif axis == 'observation':
            axis = 1
        else:
            raise UnknownAxisError(axis)

        matrix_sum = np.squeeze(np.asarray(self._data.sum(axis=axis)))

        # We only want to return a scalar if the whole matrix was summed.
        if axis is not None and matrix_sum.shape == ():
            matrix_sum = matrix_sum.reshape(1)

        return matrix_sum

    def transpose(self):
        """Transpose the contingency table

        The returned table will be an entirely new table, including copies of
        the (transposed) data, sample/observation IDs and metadata.

        Returns
        -------
        Table
            Return a new table that is the transpose of caller table.
        """
        sample_md_copy = deepcopy(self.metadata())
        obs_md_copy = deepcopy(self.metadata(axis='observation'))

        if self._data.getformat() == 'lil':
            # lil's transpose method doesn't have the copy kwarg, but all of
            # the others do.
            self._data = self._data.tocsr()

        # sample ids and observations are reversed becuase we trasposed
        return self.__class__(self._data.transpose(copy=True),
                              self.ids()[:], self.ids(axis='observation')[:],
                              sample_md_copy, obs_md_copy, self.table_id)

    def head(self, n=5, m=5):
        """Get the first n rows and m columns from self

        Parameters
        ----------
        n : int, optional
            The number of rows (observations) to get. This number must be
            greater than 0. If not specified, 5 rows will be retrieved.

        m : int, optional
            The number of columns (samples) to get. This number must be
            greater than 0. If not specified, 5 columns will be
            retrieved.

        Notes
        -----
        Like `head` for Linux like systems, requesting more rows (or columns)
        than exists will silently work.

        Raises
        ------
        IndexError
            If `n` or `m` are <= 0.

        Returns
        -------
        Table
            The subset table.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table
        >>> data = np.arange(100).reshape(5, 20)
        >>> obs_ids = ['O%d' % i for i in range(1, 6)]
        >>> samp_ids = ['S%d' % i for i in range(1, 21)]
        >>> table = Table(data, obs_ids, samp_ids)
        >>> print(table.head())  # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3  S4  S5
        O1  0.0 1.0 2.0 3.0 4.0
        O2  20.0 21.0 22.0 23.0 24.0
        O3  40.0 41.0 42.0 43.0 44.0
        O4  60.0 61.0 62.0 63.0 64.0
        O5  80.0 81.0 82.0 83.0 84.0

        """
        if n <= 0:
            raise IndexError("n cannot be <= 0.")

        if m <= 0:
            raise IndexError("m cannot be <= 0.")

        row_ids = self.ids(axis='observation')[:n]
        col_ids = self.ids(axis='sample')[:m]

        table = self.filter(row_ids, axis='observation', inplace=False)
        return table.filter(col_ids, axis='sample')

    def group_metadata(self, axis='sample'):
        """Return the group metadata of the given axis

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            Axis to search for the group metadata. Defaults to 'sample'

        Returns
        -------
        dict
            The corresponding group metadata for the given axis

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table, with group observation metadata and no group
        sample metadata:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> group_observation_md = {'tree': ('newick', '(O1:0.3,O2:0.4);')}
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               observation_group_metadata=group_observation_md)

        Get the observation group metadata:

        >>> table.group_metadata(axis='observation')
        {'tree': ('newick', '(O1:0.3,O2:0.4);')}

        Get the sample group metadata:

        >> table.group_metadata()
        None
        """
        if axis == 'sample':
            return self._sample_group_metadata
        elif axis == 'observation':
            return self._observation_group_metadata
        else:
            raise UnknownAxisError(axis)

    def ids(self, axis='sample'):
        """Return the ids along the given axis

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            Axis to return ids from. Defaults to 'sample'

        Returns
        -------
        1-D numpy array
            The ids along the given axis

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Get the ids along the observation axis:

        >>> print(table.ids(axis='observation'))
        ['O1' 'O2']

        Get the ids along the sample axis:

        >>> print(table.ids())
        ['S1' 'S2' 'S3']
        """
        if axis == 'sample':
            return self._sample_ids
        elif axis == 'observation':
            return self._observation_ids
        else:
            raise UnknownAxisError(axis)

    def update_ids(self, id_map, axis='sample', strict=True, inplace=True):
        """Update the ids along the given axis.

        Parameters
        ----------
        id_map : dict
            Mapping of old to new ids. All keys and values in this dict should
            be strings.
        axis : {'sample', 'observation'}, optional
            Axis to search for `id`. Defaults to 'sample'
        strict : bool, optional
            If ``True``, raise an error if an id is present in the given axis
            but is not a key in ``id_map``. If False, retain old identifier
            for ids that are present in the given axis but are not keys in
            ``id_map``.
        inplace : bool, optional
            If ``True`` the ids are updated in ``self``; if ``False`` the ids
            are updated in a new table is returned.

        Returns
        -------
        Table
            Table object where ids have been updated.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        TableException
            If an id from ``self`` is not in ``id_map`` and ``strict`` is
            ``True``.

        Examples
        --------
        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Define a mapping of old to new sample ids:

        >>> id_map = {'S1':'s1.1', 'S2':'s2.2', 'S3':'s3.3'}

        Get the ids along the sample axis in the table:

        >>> print(table.ids(axis='sample'))
        ['S1' 'S2' 'S3']

        Update the sample ids and get the ids along the sample axis in the
        updated table:

        >>> updated_table = table.update_ids(id_map, axis='sample')
        >>> print(updated_table.ids(axis='sample'))
        ['s1.1' 's2.2' 's3.3']
        """
        max_str_len = max([len(v) for v in id_map.values()])
        if not strict:
            ids = self.ids(axis=axis)
            max_str_len = max(max_str_len, max([len(i) for i in ids]))

        str_dtype = 'U%d' % max_str_len
        updated_ids = zeros(self.ids(axis=axis).size, dtype=str_dtype)
        for idx, old_id in enumerate(self.ids(axis=axis)):
            if strict and old_id not in id_map:
                raise TableException(
                    "Mapping not provided for %s identifier: %s. If this "
                    "identifier should not be updated, pass strict=False."
                    % (axis, old_id))

            updated_ids[idx] = id_map.get(old_id, old_id)

        # see issue #892, this protects against modifying inplace with bad
        # duplicate identifiers
        if inplace:
            if len(updated_ids) != len(set(updated_ids)):
                raise TableException("Duplicate IDs observed")

        # prepare the result object and update the ids along the specified
        # axis
        result = self if inplace else self.copy()
        if axis == 'sample':
            result._sample_ids = updated_ids
        else:
            result._observation_ids = updated_ids

        result._index_ids(None, None)

        # check for errors (specifically, we want to ensure that duplicate
        # ids haven't been introduced)
        errcheck(result)

        return result

    def _get_sparse_data(self, axis='sample'):
        """Returns the internal data in the correct sparse representation

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            Axis to search for `id`. Defaults to 'sample'

        Returns
        -------
        sparse matrix
            The data in csc (axis='sample') or csr (axis='observation')
            representation
        """
        if axis == 'sample':
            return self._data.tocsc()
        elif axis == 'observation':
            return self._data.tocsr()
        else:
            raise UnknownAxisError(axis)

    def metadata(self, id=None, axis='sample'):
        """Return the metadata of the identified sample/observation.

        Parameters
        ----------
        id : str
            ID of the sample or observation whose index will be returned.
        axis : {'sample', 'observation'}
            Axis to search for `id`.

        Returns
        -------
        defaultdict or None
            The corresponding metadata ``defaultdict`` or ``None`` of that axis
            does not have metadata.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        UnknownIDError
            If provided an unrecognized sample/observation ID.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table, with observation metadata and no sample
        metadata:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               [{'foo': 'bar'}, {'x': 'y'}], None)

        Get the metadata of the observation with ID "O2":

        >>> # casting to `dict` as the return is `defaultdict`
        >>> dict(table.metadata('O2', 'observation'))
        {'x': 'y'}

        Get the metadata of the sample with ID "S1":

        >>> table.metadata('S1', 'sample') is None
        True
        """
        if axis == 'sample':
            md = self._sample_metadata
        elif axis == 'observation':
            md = self._observation_metadata
        else:
            raise UnknownAxisError(axis)

        if id is None:
            return md

        idx = self.index(id, axis=axis)

        return md[idx] if md is not None else None

    def index(self, id, axis):
        """Return the index of the identified sample/observation.

        Parameters
        ----------
        id : str
            ID of the sample or observation whose index will be returned.
        axis : {'sample', 'observation'}
            Axis to search for `id`.

        Returns
        -------
        int
            Index of the sample/observation identified by `id`.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        UnknownIDError
            If provided an unrecognized sample/observation ID.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Get the index of the observation with ID "O2":

        >>> table.index('O2', 'observation')
        1

        Get the index of the sample with ID "S1":

        >>> table.index('S1', 'sample')
        0
        """
        idx_lookup = self._index(axis=axis)

        if id not in idx_lookup:
            raise UnknownIDError(id, axis)

        return idx_lookup[id]

    def get_value_by_ids(self, obs_id, samp_id):
        """Return value in the matrix corresponding to ``(obs_id, samp_id)``

        Parameters
        ----------
        obs_id : str
            The ID of the observation
        samp_id : str
            The ID of the sample

        Returns
        -------
        float
            The data value corresponding to the specified matrix position

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'Z3'])

        Retrieve the number of counts for observation `O1` in sample `Z3`.

        >>> print(table.get_value_by_ids('O2', 'Z3'))
        42.0

        See Also
        --------
        Table.data
        """
        return self[self.index(obs_id, 'observation'),
                    self.index(samp_id, 'sample')]

    def __str__(self):
        """Stringify self

        Default str output for a Table is just row/col ids and data values
        """
        return self.delimited_self()

    def __repr__(self):
        """Returns a high-level summary of the table's properties

        Returns
        -------
        str
            A string detailing the shape, class, number of nonzero entries, and
            table density
        """
        rows, cols = self.shape
        return '%d x %d %s with %d nonzero entries (%d%% dense)' % (
            rows, cols, repr(self.__class__), self.nnz,
            self.get_table_density() * 100
        )

    def exists(self, id, axis="sample"):
        """Returns whether id exists in axis

        Parameters
        ----------
        id: str
            id to check if exists
        axis : {'sample', 'observation'}, optional
            The axis to check

        Returns
        -------
        bool
            ``True`` if `id` exists, ``False`` otherwise

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Check whether sample ID is in the table:

        >>> table.exists('S1')
        True
        >>> table.exists('S4')
        False

        Check whether an observation ID is in the table:

        >>> table.exists('O1', 'observation')
        True
        >>> table.exists('O3', 'observation')
        False
        """
        return id in self._index(axis=axis)

    def delimited_self(self, delim='\t', header_key=None, header_value=None,
                       metadata_formatter=str,
                       observation_column_name='#OTU ID', direct_io=None):
        """Return self as a string in a delimited form

        Default str output for the Table is just row/col ids and table data
        without any metadata

        Including observation metadata in output: If ``header_key`` is not
        ``None``, the observation metadata with that name will be included
        in the delimited output. If ``header_value`` is also not ``None``, the
        observation metadata will use the provided ``header_value`` as the
        observation metadata name (i.e., the column header) in the delimited
        output.

        ``metadata_formatter``: a function which takes a metadata entry and
        returns a formatted version that should be written to file

        ``observation_column_name``: the name of the first column in the output
        table, corresponding to the observation IDs. For example, the default
        will look something like:

            #OTU ID\tSample1\tSample2
            OTU1\t10\t2
            OTU2\t4\t8
        """
        def to_utf8(i):
            if isinstance(i, bytes):
                return i.decode('utf8')
            else:
                return str(i)

        if self.is_empty():
            raise TableException("Cannot delimit self if I don't have data...")

        samp_ids = delim.join([to_utf8(i) for i in self.ids()])
        # 17 hrs of straight programming later...
        if header_key is not None:
            if header_value is None:
                raise TableException(
                    "You need to specify both header_key and header_value")

        if header_value is not None:
            if header_key is None:
                raise TableException(
                    "You need to specify both header_key and header_value")

        if header_value:
            output = [
                '# Constructed from biom file',
                f'{observation_column_name}{delim}{samp_ids}\t{header_value}'
            ]
        else:
            output = ['# Constructed from biom file',
                      f'{observation_column_name}{delim}{samp_ids}']

        if direct_io is not None:
            direct_io.writelines([i+"\n" for i in output])

        obs_metadata = self.metadata(axis='observation')
        iterable = self.ids(axis='observation')
        end_line = '' if direct_io is None else '\n'

        for obs_id, obs_values in zip(iterable,
                                      self._iter_obs()):
            str_obs_vals = delim.join(map(str, self._to_dense(obs_values)))
            obs_id = to_utf8(obs_id)
            if header_key and obs_metadata is not None:
                md = obs_metadata[self._obs_index[obs_id]]
                md_out = metadata_formatter(md.get(header_key, None))
                output_row = '%s%s%s\t%s%s' % \
                    (obs_id, delim, str_obs_vals, md_out, end_line)

                if direct_io is None:
                    output.append(output_row)
                else:
                    direct_io.write(output_row)
            else:
                output_row = '%s%s%s%s' % \
                            (obs_id, delim, str_obs_vals, end_line)
                if direct_io is None:
                    output.append(output_row)
                else:
                    direct_io.write(output_row)

        return '\n'.join(output)

    def is_empty(self):
        """Check whether the table is empty

        Returns
        -------
        bool
            ``True`` if the table is empty, ``False`` otherwise
        """
        if not self.ids().size or not self.ids(axis='observation').size:
            return True
        else:
            return False

    def __iter__(self):
        """See ``biom.table.Table.iter``"""
        return self.iter()

    def _iter_samp(self):
        """Return sample vectors of data matrix vectors"""
        for c in range(self.shape[1]):
            # this pulls out col vectors but need to convert to the expected
            # row vector
            colvec = self._get_col(c)
            yield colvec.transpose(copy=True)

    def _iter_obs(self):
        """Return observation vectors of data matrix"""
        for r in range(self.shape[0]):
            yield self._get_row(r)

    def get_table_density(self):
        """Returns the fraction of nonzero elements in the table.

        Returns
        -------
        float
            The fraction of nonzero elements in the table
        """
        density = 0.0

        if not self.is_empty():
            density = (self.nnz /
                       (len(self.ids()) * len(self.ids(axis='observation'))))

        return density

    def descriptive_equality(self, other):
        """For use in testing, describe how the tables are not equal"""
        if not isinstance(other, self.__class__):
            return "Tables are not of comparable classes"
        if not self.type == other.type:
            return "Tables are not the same type"
        if not np.array_equal(self.ids(axis='observation'),
                              other.ids(axis='observation')):
            return "Observation IDs are not the same"
        if not np.array_equal(self.ids(), other.ids()):
            return "Sample IDs are not the same"
        if not np.array_equal(self.metadata(axis='observation'),
                              other.metadata(axis='observation')):
            return "Observation metadata are not the same"
        if not np.array_equal(self.metadata(), other.metadata()):
            return "Sample metadata are not the same"
        if not self._data_equality(other._data):
            return "Data elements are not the same"

        return "Tables appear equal"

    def __eq__(self, other):
        """Equality is determined by the data matrix, metadata, and IDs"""
        if not isinstance(other, self.__class__):
            return False
        if self.type != other.type:
            return False
        if not np.array_equal(self.ids(axis='observation'),
                              other.ids(axis='observation')):
            return False
        if not np.array_equal(self.ids(), other.ids()):
            return False
        if not np.array_equal(self.metadata(axis='observation'),
                              other.metadata(axis='observation')):
            return False
        if not np.array_equal(self.metadata(), other.metadata()):
            return False
        if not self._data_equality(other._data):
            return False

        return True

    def _data_equality(self, other):
        """Return ``True`` if both matrices are equal.

        Matrices are equal iff the following items are equal:
        - shape
        - dtype
        - size (nnz)
        - matrix data (more expensive, so checked last)

        The sparse format does not need to be the same between the two
        matrices. ``self`` and ``other`` will be converted to csr format if
        necessary before performing the final comparison.

        """
        if self._data.shape != other.shape:
            return False

        if self._data.dtype != other.dtype:
            return False

        if self._data.nnz != other.nnz:
            return False

        self._data = self._data.tocsr()
        other = other.tocsr()

        if (self._data != other).nnz > 0:
            return False

        return True

    def __ne__(self, other):
        return not (self == other)

    def data(self, id, axis='sample', dense=True):
        """Returns data associated with an `id`

        Parameters
        ----------
        id : str
            ID of the sample or observation whose data will be returned.
        axis : {'sample', 'observation'}
            Axis to search for `id`.
        dense : bool, optional
            If ``True``, return data as dense

        Returns
        -------
        np.ndarray or scipy.sparse.spmatrix
            np.ndarray if ``dense``, otherwise scipy.sparse.spmatrix

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> from biom import example_table
        >>> example_table.data('S1', axis='sample')
        array([ 0.,  3.])

        See Also
        --------
        Table.get_value_by_ids

        """
        if axis == 'sample':
            data = self[:, self.index(id, 'sample')]
        elif axis == 'observation':
            data = self[self.index(id, 'observation'), :]
        else:
            raise UnknownAxisError(axis)

        if dense:
            return self._to_dense(data)
        else:
            return data

    def copy(self):
        """Returns a copy of the table"""
        return self.__class__(self._data.copy(),
                              self.ids(axis='observation').copy(),
                              self.ids().copy(),
                              deepcopy(self.metadata(axis='observation')),
                              deepcopy(self.metadata()),
                              self.table_id,
                              type=self.type)

    def iter_data(self, dense=True, axis='sample'):
        """Yields axis values

        Parameters
        ----------
        dense : bool, optional
            Defaults to ``True``. If ``False``, yield compressed sparse row or
            compressed sparse columns if `axis` is 'observation' or 'sample',
            respectively.
        axis : {'sample', 'observation'}, optional
            Axis to iterate over.

        Returns
        -------
        generator
            Yields list of values for each value in `axis`

        Raises
        ------
        UnknownAxisError
            If axis other than 'sample' or 'observation' passed

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table
        >>> data = np.arange(30).reshape(3,10) # 3 X 10 OTU X Sample table
        >>> obs_ids = ['o1', 'o2', 'o3']
        >>> sam_ids = ['s%i' %i for i in range(1,11)]
        >>> bt = Table(data, observation_ids=obs_ids, sample_ids=sam_ids)

        Lets find the sample with the largest sum

        >>> sample_gen = bt.iter_data(axis='sample')
        >>> max_sample_count = max([sample.sum() for sample in sample_gen])
        >>> print(max_sample_count)
        57.0
        """
        if axis == "sample":
            for samp_v in self._iter_samp():
                if dense:
                    yield self._to_dense(samp_v)
                else:
                    yield samp_v
        elif axis == "observation":
            for obs_v in self._iter_obs():
                if dense:
                    yield self._to_dense(obs_v)
                else:
                    yield obs_v
        else:
            raise UnknownAxisError(axis)

    def iter(self, dense=True, axis='sample'):
        """Yields ``(value, id, metadata)``


        Parameters
        ----------
        dense : bool, optional
            Defaults to ``True``. If ``False``, yield compressed sparse row or
            compressed sparse columns if `axis` is 'observation' or 'sample',
            respectively.
        axis : {'sample', 'observation'}, optional
            The axis to iterate over.

        Returns
        -------
        GeneratorType
            A generator that yields (values, id, metadata)

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'Z3'])

        Iter over samples and keep those that start with an Z:

        >>> [(values, id, metadata)
        ...     for values, id, metadata in table.iter() if id[0]=='Z']
        [(array([  1.,  42.]), 'Z3', None)]

        Iter over observations and add the 2nd column of the values

        >>> col = [values[1] for values, id, metadata in table.iter()]
        >>> sum(col)
        46.0
        """
        ids = self.ids(axis=axis)
        metadata = self.metadata(axis=axis)
        if axis == 'sample':
            iter_ = self._iter_samp()
        elif axis == 'observation':
            iter_ = self._iter_obs()
        else:
            raise UnknownAxisError(axis)

        if metadata is None:
            metadata = (None,) * len(ids)

        iter_ = self.iter_data(axis=axis, dense=dense)

        return zip(iter_, ids, metadata)

    def iter_pairwise(self, dense=True, axis='sample', tri=True, diag=False):
        """Pairwise iteration over self

        Parameters
        ----------
        dense : bool, optional
            Defaults to ``True``. If ``False``, yield compressed sparse row or
            compressed sparse columns if `axis` is 'observation' or 'sample',
            respectively.
        axis : {'sample', 'observation'}, optional
            The axis to iterate over.
        tri : bool, optional
            If ``True``, just yield [i, j] and not [j, i]
        diag : bool, optional
            If ``True``, yield [i, i]

        Returns
        -------
        GeneratorType
            Yields [(val_i, id_i, metadata_i), (val_j, id_j, metadata_j)]

        Raises
        ------
        UnknownAxisError

        Examples
        --------
        >>> from biom import example_table

        By default, only the upper triangle without the diagonal  of the
        resulting pairwise combinations is yielded.

        >>> iter_ = example_table.iter_pairwise()
        >>> for (val_i, id_i, md_i), (val_j, id_j, md_j) in iter_:
        ...     print(id_i, id_j)
        S1 S2
        S1 S3
        S2 S3

        The full pairwise combinations can also be yielded though.

        >>> iter_ = example_table.iter_pairwise(tri=False, diag=True)
        >>> for (val_i, id_i, md_i), (val_j, id_j, md_j) in iter_:
        ...     print(id_i, id_j)
        S1 S1
        S1 S2
        S1 S3
        S2 S1
        S2 S2
        S2 S3
        S3 S1
        S3 S2
        S3 S3

        """
        metadata = self.metadata(axis=axis)
        ids = self.ids(axis=axis)

        if metadata is None:
            metadata = (None,) * len(ids)

        ind = np.arange(len(ids))
        diag_v = 1 - diag  # for offseting tri_f, where a 0 includes the diag

        if tri:
            def tri_f(idx):
                return ind[idx+diag_v:]
        else:
            def tri_f(idx):
                return np.hstack([ind[:idx], ind[idx+diag_v:]])

        for idx, i in enumerate(ind):
            id_i = ids[i]
            md_i = metadata[i]
            data_i = self.data(id_i, axis=axis, dense=dense)

            for j in tri_f(idx):
                id_j = ids[j]
                md_j = metadata[j]
                data_j = self.data(id_j, axis=axis, dense=dense)

                yield ((data_i, id_i, md_i), (data_j, id_j, md_j))

    def sort_order(self, order, axis='sample'):
        """Return a new table with `axis` in `order`

        Parameters
        ----------
        order : iterable
            The desired order for axis
        axis : {'sample', 'observation'}, optional
            The axis to operate on

        Returns
        -------
        Table
            A table where the observations or samples are sorted according to
            `order`

        Examples
        --------

        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[1, 0, 4], [1, 3, 0]])
        >>> table = Table(data, ['O2', 'O1'], ['S2', 'S1', 'S3'])
        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S2  S1  S3
        O2  1.0 0.0 4.0
        O1  1.0 3.0 0.0

        Sort the table using a list of samples:

        >>> sorted_table = table.sort_order(['S2', 'S3', 'S1'])
        >>> print(sorted_table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	S2	S3	S1
        O2	1.0	4.0	0.0
        O1	1.0	0.0	3.0


        Additionally you could sort the table's observations:

        >>> sorted_table = table.sort_order(['O1', 'O2'], axis="observation")
        >>> print(sorted_table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	S2	S1	S3
        O1	1.0	3.0	0.0
        O2	1.0	0.0	4.0

        """
        fancy = np.array([self.index(i, axis=axis) for i in order], dtype=int)
        metadata = self.metadata(axis=axis)
        if metadata is not None:
            metadata = np.array(metadata)[fancy]

        if axis == 'sample':
            mat = self.matrix_data[:, fancy]
            return self.__class__(mat,
                                  self.ids(axis='observation')[:], order[:],
                                  self.metadata(axis='observation'), metadata,
                                  self.table_id, self.type)

        elif axis == 'observation':
            mat = self.matrix_data[fancy, :]
            return self.__class__(mat,
                                  order[:], self.ids()[:],
                                  metadata, self.metadata(), self.table_id,
                                  self.type)
        else:
            raise UnknownAxisError(axis)

    def sort(self, sort_f=natsort, axis='sample'):
        """Return a table sorted along axis

        Parameters
        ----------
        sort_f : function, optional
            Defaults to ``biom.util.natsort``. A function that takes a list of
            values and sorts it
        axis : {'sample', 'observation'}, optional
            The axis to operate on

        Returns
        -------
        biom.Table
            A table whose samples or observations are sorted according to the
            `sort_f` function

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[1, 0, 4], [1, 3, 0]])
        >>> table = Table(data, ['O2', 'O1'], ['S2', 'S1', 'S3'])
        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S2  S1  S3
        O2  1.0 0.0 4.0
        O1  1.0 3.0 0.0

        Sort the order of samples in the table using the default natural
        sorting:

        >>> new_table = table.sort()
        >>> print(new_table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3
        O2  0.0 1.0 4.0
        O1  3.0 1.0 0.0

        Sort the order of observations in the table using the default natural
        sorting:

        >>> new_table = table.sort(axis='observation')
        >>> print(new_table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S2  S1  S3
        O1  1.0 3.0 0.0
        O2  1.0 0.0 4.0

        Sort the samples in reverse order using a custom sort function:

        >>> sort_f = lambda x: list(sorted(x, reverse=True))
        >>> new_table = table.sort(sort_f=sort_f)
        >>> print(new_table)  # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S3  S2  S1
        O2  4.0 1.0 0.0
        O1  0.0 1.0 3.0
        """
        return self.sort_order(sort_f(self.ids(axis=axis)), axis=axis)

    def filter(self, ids_to_keep, axis='sample', invert=False, inplace=True):
        """Filter a table based on a function or iterable.

        Parameters
        ----------
        ids_to_keep : iterable, or function(values, id, metadata) -> bool
            If a function, it will be called with the values of the
            sample/observation, its id (a string) and the dictionary
            of metadata of each sample/observation, and must return a
            boolean. If it's an iterable, it must be a list of ids to
            keep.
        axis : {'sample', 'observation'}, optional
            It controls whether to filter samples or observations and
            defaults to "sample".
        invert : bool, optional
            Defaults to ``False``. If set to ``True``, discard samples or
            observations where `ids_to_keep` returns True
        inplace : bool, optional
            Defaults to ``True``. Whether to return a new table or modify
            itself.

        Returns
        -------
        biom.Table
            Returns itself if `inplace`, else returns a new filtered table.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table, with observation metadata and sample
        metadata:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               [{'full_genome_available': True},
        ...                {'full_genome_available': False}],
        ...               [{'sample_type': 'a'}, {'sample_type': 'a'},
        ...                {'sample_type': 'b'}])

        Define a function to keep only samples with sample_type == 'a'. This
        will drop sample S3, which has sample_type 'b':

        >>> filter_fn = lambda val, id_, md: md['sample_type'] == 'a'

        Get a filtered version of the table, leaving the original table
        untouched:

        >>> new_table = table.filter(filter_fn, inplace=False)
        >>> print(table.ids())
        ['S1' 'S2' 'S3']
        >>> print(new_table.ids())
        ['S1' 'S2']

        Using the same filtering function, discard all samples with sample_type
        'a'. This will keep only sample S3, which has sample_type 'b':

        >>> new_table = table.filter(filter_fn, inplace=False, invert=True)
        >>> print(table.ids())
        ['S1' 'S2' 'S3']
        >>> print(new_table.ids())
        ['S3']

        Filter the table in-place using the same function (drop all samples
        where sample_type is not 'a'):

        >>> table.filter(filter_fn)
        2 x 2 <class 'biom.table.Table'> with 2 nonzero entries (50% dense)
        >>> print(table.ids())
        ['S1' 'S2']

        Filter out all observations in the table that do not have
        full_genome_available == True. This will filter out observation O2:

        >>> filter_fn = lambda val, id_, md: md['full_genome_available']
        >>> table.filter(filter_fn, axis='observation')
        1 x 2 <class 'biom.table.Table'> with 0 nonzero entries (0% dense)
        >>> print(table.ids(axis='observation'))
        ['O1']

        """
        table = self if inplace else self.copy()

        metadata = table.metadata(axis=axis)
        ids = table.ids(axis=axis)
        index = self._index(axis=axis)
        axis = table._axis_to_num(axis=axis)
        arr = table._data
        arr, ids, metadata = _filter(arr,
                                     ids,
                                     metadata,
                                     index,
                                     ids_to_keep,
                                     axis,
                                     invert=invert)

        table._data = arr
        if axis == 1:
            table._sample_ids = ids
            table._sample_metadata = metadata
            table._index_ids(self._obs_index.copy(), None)
        elif axis == 0:
            table._observation_ids = ids
            table._observation_metadata = metadata
            table._index_ids(None, self._sample_index.copy())

        errcheck(table)

        return table

    def partition(self, f, axis='sample', remove_empty=False,
                  ignore_none=False):
        """Yields partitions

        Parameters
        ----------
        f : function, dict
            `f` is given the ID and metadata of the vector and must return
            what partition the vector is part of. If `dict`, a mapping of
            either ID -> group, or group -> [list, of, ID] must be provided.
        axis : {'sample', 'observation'}, optional
            The axis to iterate over
        remove_empty : bool, optional
            If `True`, remove empty vectors from a partition. Default is
            `False`.
        ignore_none : bool, optional
            If `True`, ignore partitions with the label `None`. Default is
            `False`.

        Returns
        -------
        GeneratorType
            A generator that yields (partition, `Table`)

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table
        >>> from biom.util import unzip

        Create a 2x3 BIOM table, with observation metadata and sample
        metadata:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               [{'full_genome_available': True},
        ...                {'full_genome_available': False}],
        ...               [{'sample_type': 'a'}, {'sample_type': 'a'},
        ...                {'sample_type': 'b'}])

        Define a function to bin by sample_type

        >>> f = lambda id_, md: md['sample_type']

        Partition the table and view results

        >>> bins, tables = table.partition(f)
        >>> print(bins[1]) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2
        O1  0.0 0.0
        O2  1.0 3.0
        >>> print(tables[1]) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S3
        O1  1.0
        O2  42.0
        """
        # we are not checking for whether the IDs are or are not present as
        # that introduces complexity of `strict`. Deferring that for now.
        if isinstance(f, dict):
            test = list(f.values())[0]

            if isinstance(test, (list, tuple)):
                # group -> [list, of, ids]
                mapping = {}
                for grp, ids in f.items():
                    for id_ in ids:
                        mapping[id_] = grp

            elif isinstance(test, str):
                # id_ -> grp
                mapping = f

            else:
                raise ValueError(f"Unable to handle a type of `{type(test)}` "
                                 "with mapping")

            def part_f(i, m):
                return mapping.get(i)
        else:
            part_f = f

        partitions = {}
        # conversion of vector types is not necessary, vectors are not
        # being passed to an arbitrary function
        for vals, id_, md in self.iter(dense=False, axis=axis):
            part = part_f(id_, md)

            if ignore_none and part is None:
                continue

            # try to make it hashable...
            if not isinstance(part, Hashable):
                part = tuple(part)

            if part not in partitions:
                partitions[part] = [[], [], []]

            partitions[part][0].append(id_)
            partitions[part][1].append(vals)
            partitions[part][2].append(md)

        md = self.metadata(axis=self._invert_axis(axis))

        for part, (ids, values, metadata) in partitions.items():
            if axis == 'sample':
                data = self._conv_to_self_type(values, transpose=True)
                samp_ids = ids
                samp_md = metadata
                obs_ids = self.ids(axis='observation')[:]
                obs_md = md[:] if md is not None else None
                indices = {'observation_index': self._obs_index.copy()}

            elif axis == 'observation':
                data = self._conv_to_self_type(values, transpose=False)
                obs_ids = ids
                obs_md = metadata
                samp_ids = self.ids()[:]
                samp_md = md[:] if md is not None else None
                indices = {'sample_index': self._sample_index.copy()}

            tab = Table(data, obs_ids, samp_ids, obs_md, samp_md,
                        self.table_id, type=self.type, validate=False,
                        **indices)

            if remove_empty:
                tab.remove_empty(inplace=True)

            yield part, tab

    def collapse(self, f, collapse_f=None, norm=True, min_group_size=1,
                 include_collapsed_metadata=True, one_to_many=False,
                 one_to_many_mode='add', one_to_many_md_key='Path',
                 strict=False, axis='sample'):
        """Collapse partitions in a table by metadata or by IDs

        Partition data by metadata or IDs and then collapse each partition into
        a single vector.

        If `include_collapsed_metadata` is ``True``, the metadata for the
        collapsed partition will be a category named 'collapsed_ids', in which
        a list of the original ids that made up the partition is retained

        The remainder is only relevant to setting `one_to_many` to ``True``.

        If `one_to_many` is ``True``, allow vectors to collapse into multiple
        bins if the metadata describe a one-many relationship. Supplied
        functions must allow for iteration support over the metadata key and
        must return a tuple of (path, bin) as to describe both the path in the
        hierarchy represented and the specific bin being collapsed into. The
        uniqueness of the bin is _not_ based on the path but by the name of the
        bin.

        The metadata value for the corresponding collapsed column may include
        more (or less) information about the collapsed data. For example, if
        collapsing "FOO", and there are vectors that span three associations A,
        B, and C, such that vector 1 spans A and B, vector 2 spans B and C and
        vector 3 spans A and C, the resulting table will contain three
        collapsed vectors:

        - A, containing original vectors 1 and 3
        - B, containing original vectors 1 and 2
        - C, containing original vectors 2 and 3

        If a vector maps to the same partition multiple times, it will be
        counted multiple times.

        There are two supported modes for handling one-to-many relationships
        via `one_to_many_mode`: ``add`` and `divide`. ``add`` will add the
        vector counts to each partition that the vector maps to, which may
        increase the total number of counts in the output table. ``divide``
        will divide a vectors's counts by the number of metadata that the
        vector has before adding the counts to each partition. This will not
        increase the total number of counts in the output table.

        If `one_to_many_md_key` is specified, that becomes the metadata
        key that describes the collapsed path. If a value is not specified,
        then it defaults to 'Path'.

        If `strict` is specified, then all metadata pathways operated on
        must be indexable by `metadata_f`.

        `one_to_many` and `norm` are not supported together.

        `one_to_many` and `collapse_f` are not supported together.

        `one_to_many` and `min_group_size` are not supported together.

        A final note on space consumption. At present, the `one_to_many`
        functionality requires a temporary dense matrix representation.

        Parameters
        ----------
        f : function
            Function that is used to determine what partition a vector belongs
            to
        collapse_f : function, optional
            Function that collapses a partition in a one-to-one collapse. The
            expected function signature is:

                dense or sparse_vector <- collapse_f(Table, axis)

            Defaults to a pairwise add.

        norm : bool, optional
            Defaults to ``True``. If ``True``, normalize the resulting table
        min_group_size : int, optional
            Defaults to ``1``. The minimum size of a partition when performing
            a one-to-one collapse
        include_collapsed_metadata : bool, optional
            Defaults to ``True``. If ``True``, retain the collapsed metadata
            keyed by the original IDs of the associated vectors
        one_to_many : bool, optional
            Defaults to ``False``. Perform a one-to-many collapse
        one_to_many_mode : {'add', 'divide'}, optional
            The way to reduce two vectors in a one-to-many collapse
        one_to_many_md_key : str, optional
            Defaults to "Path". If `include_collapsed_metadata` is ``True``,
            store the original vector metadata under this key
        strict : bool, optional
            Defaults to ``False``. Requires full pathway data within a
            one-to-many structure
        axis : {'sample', 'observation'}, optional
            The axis to collapse

        Returns
        -------
        Table
            The collapsed table

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a ``Table``

        >>> dt_rich = Table(
        ...    np.array([[5, 6, 7], [8, 9, 10], [11, 12, 13]]),
        ...    ['1', '2', '3'], ['a', 'b', 'c'],
        ...    [{'taxonomy': ['k__a', 'p__b']},
        ...     {'taxonomy': ['k__a', 'p__c']},
        ...     {'taxonomy': ['k__a', 'p__c']}],
        ...    [{'barcode': 'aatt'},
        ...     {'barcode': 'ttgg'},
        ...     {'barcode': 'aatt'}])
        >>> print(dt_rich) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID a   b   c
        1   5.0 6.0 7.0
        2   8.0 9.0 10.0
        3   11.0    12.0    13.0

        Create Function to determine what partition a vector belongs to

        >>> bin_f = lambda id_, x: x['taxonomy'][1]
        >>> obs_phy = dt_rich.collapse(
        ...    bin_f, norm=False, min_group_size=1,
        ...    axis='observation').sort(axis='observation')
        >>> print(obs_phy) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID a   b   c
        p__b    5.0 6.0 7.0
        p__c    19.0    21.0    23.0
        """
        collapsed_data = []
        collapsed_ids = []

        if include_collapsed_metadata:
            collapsed_md = []
        else:
            collapsed_md = None

        if one_to_many_mode not in ['add', 'divide']:
            raise ValueError("Unrecognized one-to-many mode '%s'. Must be "
                             "either 'add' or 'divide'." % one_to_many_mode)

        # transpose is only necessary in the one-to-one case
        # new_data_shape is only necessary in the one-to-many case
        # axis_update is only necessary in the one-to-many case
        def axis_ids_md(t):
            return (t.ids(axis=axis), t.metadata(axis=axis))

        if axis == 'sample':
            transpose = True

            def axis_update(offaxis, onaxis):
                return (offaxis, onaxis)

        elif axis == 'observation':
            transpose = False

            def axis_update(offaxis, onaxis):
                return (onaxis, offaxis)

        else:
            raise UnknownAxisError(axis)

        if one_to_many:
            if norm:
                raise AttributeError(
                    "norm and one_to_many are not supported together")

            # determine the collapsed pathway
            # we drop all other associated metadata
            new_md = {}
            md_count = {}

            for id_, md in zip(*axis_ids_md(self)):
                md_iter = f(id_, md)
                num_md = 0
                while True:
                    try:
                        pathway, partition = next(md_iter)
                    except IndexError:
                        # if a pathway is incomplete
                        if strict:
                            # bail if strict
                            err = "Incomplete pathway, ID: %s, metadata: %s" %\
                                  (id_, md)
                            raise IndexError(err)
                        else:
                            # otherwise ignore
                            continue
                    except StopIteration:
                        break

                    new_md[partition] = pathway
                    num_md += 1

                md_count[id_] = num_md

            idx_lookup = {part: i for i, part in enumerate(sorted(new_md))}

            # We need to store floats, not ints, as things won't always divide
            # evenly.
            dtype = np.float64 if one_to_many_mode == 'divide' else self.dtype

            if axis == 'observation':
                new_data = dok_matrix((len(self.ids(axis='sample')),
                                       len(new_md)),
                                      dtype=dtype)
            else:
                new_data = dok_matrix((len(self.ids(axis='observation')),
                                       len(new_md)), dtype=dtype)

            # for each vector
            # for each bin in the metadata
            # for each partition associated with the vector
            for vals, id_, md in self.iter(axis=axis, dense=False):
                md_iter = f(id_, md)

                while True:
                    try:
                        pathway, part = next(md_iter)
                    except IndexError:
                        # if a pathway is incomplete
                        if strict:
                            # bail if strict, should never get here...
                            err = "Incomplete pathway, ID: %s, metadata: %s" %\
                                  (id_, md)
                            raise IndexError(err)
                        else:
                            # otherwise ignore
                            continue
                    except StopIteration:
                        break

                    # TODO: refactor Table.collapse(..., one_to_many=True) so
                    # writes into new_data are performed without regard to
                    # the requested axis, and perform a single transpose at the
                    # end. Right now we incur many calls to `axis_update` which
                    # could be avoided. However, this refactor is likely
                    # complex to do correctly, so punting for now as we don't
                    # yet have data showing this is a real world performance
                    # concern.
                    column = idx_lookup[part]
                    if one_to_many_mode == 'add':
                        for vidx, v in zip(vals.indices, vals.data):
                            new_data[vidx, column] += v
                    else:
                        dv = md_count[id_]
                        tmp = vals / dv
                        for vidx, v in zip(tmp.indices, tmp.data):
                            new_data[vidx, column] += v

            if include_collapsed_metadata:
                # reassociate pathway information
                for k, i in sorted(idx_lookup.items(), key=itemgetter(1)):
                    collapsed_md.append({one_to_many_md_key: new_md[k]})

            # get the new sample IDs
            collapsed_ids = [k for k, i in sorted(idx_lookup.items(),
                                                  key=itemgetter(1))]

            # convert back to self type
            if axis == 'observation':
                new_data = csr_matrix(new_data.T)
            else:
                new_data = csc_matrix(new_data)

            data = self._conv_to_self_type(new_data)
        else:
            if collapse_f is None:
                def collapse_f(t, axis):
                    return t.sum(axis)

            for part, table in self.partition(f, axis=axis):
                axis_ids, axis_md = axis_ids_md(table)

                if len(axis_ids) < min_group_size:
                    continue

                redux_data = collapse_f(table, self._invert_axis(axis))
                if norm:
                    redux_data /= len(axis_ids)

                collapsed_data.append(self._conv_to_self_type(redux_data))
                collapsed_ids.append(part)

                if include_collapsed_metadata:
                    # retain metadata but store by original id
                    collapsed_md.append({'collapsed_ids': axis_ids.tolist()})

            data = self._conv_to_self_type(collapsed_data, transpose=transpose)

        # if the table is empty
        errcheck(self, 'empty')

        md = self.metadata(axis=self._invert_axis(axis))
        if axis == 'sample':
            sample_ids = collapsed_ids
            sample_md = collapsed_md
            obs_ids = self.ids(axis='observation')[:]
            obs_md = md if md is not None else None
        else:
            sample_ids = self.ids()[:]
            obs_ids = collapsed_ids
            obs_md = collapsed_md
            sample_md = md if md is not None else None

        return Table(data, obs_ids, sample_ids, obs_md, sample_md,
                     self.table_id, type=self.type)

    def _invert_axis(self, axis):
        """Invert an axis"""
        if axis == 'sample':
            return 'observation'
        elif axis == 'observation':
            return 'sample'
        else:
            return UnknownAxisError(axis)

    def _axis_to_num(self, axis):
        """Convert str axis to numerical axis"""
        if axis == 'sample':
            return 1
        elif axis == 'observation':
            return 0
        else:
            raise UnknownAxisError(axis)

    def min(self, axis='sample'):
        """Get the minimum nonzero value over an axis

        Parameters
        ----------
        axis : {'sample', 'observation', 'whole'}, optional
            Defaults to "sample". The axis over which to calculate minima.

        Returns
        -------
        scalar of self.dtype or np.array of self.dtype

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> from biom import example_table
        >>> print(example_table.min(axis='sample'))
        [ 3.  1.  2.]

        """
        if axis not in ['sample', 'observation', 'whole']:
            raise UnknownAxisError(axis)

        if axis == 'whole':
            min_val = np.inf
            for data in self.iter_data(dense=False):
                # only min over the actual nonzero values
                min_val = min(min_val, data.data.min())
        else:
            min_val = zeros(len(self.ids(axis=axis)), dtype=self.dtype)

            for idx, data in enumerate(self.iter_data(dense=False, axis=axis)):
                min_val[idx] = data.data.min()

        return min_val

    def max(self, axis='sample'):
        """Get the maximum nonzero value over an axis

        Parameters
        ----------
        axis : {'sample', 'observation', 'whole'}, optional
            Defaults to "sample". The axis over which to calculate maxima.

        Returns
        -------
        scalar of self.dtype or np.array of self.dtype

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> from biom import example_table
        >>> print(example_table.max(axis='observation'))
        [ 2.  5.]

        """
        if axis not in ['sample', 'observation', 'whole']:
            raise UnknownAxisError(axis)

        if axis == 'whole':
            max_val = -np.inf
            for data in self.iter_data(dense=False):
                # only min over the actual nonzero values
                max_val = max(max_val, data.data.max())
        else:
            max_val = np.empty(len(self.ids(axis=axis)), dtype=self.dtype)

            for idx, data in enumerate(self.iter_data(dense=False, axis=axis)):
                max_val[idx] = data.data.max()

        return max_val

    def subsample(self, n, axis='sample', by_id=False, with_replacement=False,
                  seed=None):
        """Randomly subsample without replacement.

        Parameters
        ----------
        n : int
            Number of items to subsample from `counts`.
        axis : {'sample', 'observation'}, optional
            The axis to sample over
        by_id : boolean, optional
            If `False`, the subsampling is based on the counts contained in the
            matrix (e.g., rarefaction). If `True`, the subsampling is based on
            the IDs (e.g., fetch a random subset of samples). Default is
            `False`.
        with_replacement : boolean, optional
            If `False` (default), subsample without replacement. If `True`,
            resample with replacement via the multinomial distribution.
            Should not be `True` if `by_id` is `True`. Important: If `True`,
            samples with a sum below `n` are retained.
        seed : int, optional
            If provided, set the numpy random seed with this value

        Returns
        -------
        biom.Table
            A subsampled version of self

        Raises
        ------
        ValueError
            - If `n` is less than zero.
            - If `by_id` and `with_replacement` are both True.

        Notes
        -----
        If subsampling is performed without replacement, vectors with a sum
        less than `n` are omitted from the result. This condition is not held
        when operating with replacement.

        This code assumes absolute abundance if `by_id` is False.

        If subsampling with replacement, `np.ceil` is applied prior to
        calculating p-values to ensure that low-abundance features have a
        chance to be sampled.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table
        >>> table = Table(np.array([[0, 2, 3], [1, 0, 2]]), ['O1', 'O2'],
        ...               ['S1', 'S2', 'S3'])

        Subsample 1 item over the sample axis by value (e.g., rarefaction):

        >>> print(table.subsample(1).sum(axis='sample'))
        [ 1.  1.  1.]

        Subsample 2 items over the sample axis, note that 'S1' is filtered out:

        >>> ss = table.subsample(2)
        >>> print(ss.sum(axis='sample'))
        [ 2.  2.]
        >>> print(ss.ids())
        ['S2' 'S3']

        Subsample by IDs over the sample axis. For this example, we're going to
        randomly select 2 samples and do this 100 times, and then print out the
        set of IDs observed.

        >>> ids = set([tuple(table.subsample(2, by_id=True).ids())
        ...            for i in range(100)])
        >>> print(sorted(ids))
        [('S1', 'S2'), ('S1', 'S3'), ('S2', 'S3')]

        """
        if n < 0:
            raise ValueError("n cannot be negative.")

        if with_replacement and by_id:
            raise ValueError("by_id and with_replacement cannot both be True")

        table = self.copy()

        rng = np.random.default_rng(seed)

        if by_id:
            ids = table.ids(axis=axis).copy()
            rng.shuffle(ids)
            subset = set(ids[:n])
            table.filter(lambda v, i, md: i in subset, axis=axis)
        else:
            data = table._get_sparse_data()
            subsample(data, n, with_replacement, rng)
            table._data = data

            table.filter(lambda v, i, md: v.sum() > 0, axis=axis)

        inv_axis = self._invert_axis(axis)
        table.filter(lambda v, i, md: v.sum() > 0, axis=inv_axis)

        return table

    def pa(self, inplace=True):
        """Convert the table to presence/absence data

        Parameters
        ----------
        inplace : bool, optional
            Defaults to ``True``

        Returns
        -------
        Table
            Returns itself if `inplace`, else returns a new presence/absence
            table.

        Examples
        --------
        >>> from biom.table import Table
        >>> import numpy as np

        Create a 2x3 BIOM table

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Convert to presence/absence data

        >>> _ = table.pa()
        >>> print(table.data('O1', 'observation'))
        [ 0.  0.  1.]
        >>> print(table.data('O2', 'observation'))
        [ 1.  1.  1.]
        """
        def transform_f(data, id_, metadata):
            return np.where(data != 0, 1., 0.)

        return self.transform(transform_f, inplace=inplace)

    def transform(self, f, axis='sample', inplace=True):
        """Iterate over `axis`, applying a function `f` to each vector.

        Only non null values can be modified and the density of the
        table can't increase. However, zeroing values is fine.

        Parameters
        ----------
        f : function(data, id, metadata) -> new data
            A function that takes three values: an array of nonzero
            values corresponding to each observation or sample, an
            observation or sample id, and an observation or sample
            metadata entry. It must return an array of transformed
            values that replace the original values.
        axis : {'sample', 'observation'}, optional
            The axis to operate on. Can be "sample" or "observation".
        inplace : bool, optional
            Defaults to ``True``. Whether to return a new table or modify
            itself.

        Returns
        -------
        biom.Table
            Returns itself if `inplace`, else returns a new transformed table.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 table

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               [{'foo': 'bar'}, {'x': 'y'}], None)
        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3
        O1  0.0 0.0 1.0
        O2  1.0 3.0 42.0

        Create a transform function

        >>> f = lambda data, id_, md: data / 2

        Transform to a new table on samples

        >>> table2 = table.transform(f, 'sample', False)
        >>> print(table2) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3
        O1  0.0 0.0 0.5
        O2  0.5 1.5 21.0

        `table` hasn't changed

        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3
        O1  0.0 0.0 1.0
        O2  1.0 3.0 42.0

        Tranform in place on observations

        >>> table3 = table.transform(f, 'observation', True)

        `table` is different now

        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3
        O1  0.0 0.0 0.5
        O2  0.5 1.5 21.0

        but the table returned (`table3`) is the same as `table`

        >>> print(table3) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2  S3
        O1  0.0 0.0 0.5
        O2  0.5 1.5 21.0

        """
        table = self if inplace else self.copy()

        metadata = table.metadata(axis=axis)
        ids = table.ids(axis=axis)
        arr = table._get_sparse_data(axis=axis)

        axis = table._axis_to_num(axis)

        _transform(arr, ids, metadata, f, axis)
        arr.eliminate_zeros()

        table._data = arr

        return table

    def rankdata(self, axis='sample', inplace=True, method='average'):
        """Convert values to rank abundances from smallest to largest

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            The axis to use for ranking.
        inplace : bool, optional
            Defaults to ``True``. If ``True``, performs the ranking in
            place. Otherwise, returns a new table with ranking applied.
        method : str, optional
            The method for handling ties in counts. This can be any valid
            string that can be passed to `scipy.stats.rankdata`.

        Returns
        -------
        biom.Table
            The rank-abundance-transformed table.

        Raises
        ------
        ValueError
            If unknown ``method`` is provided.

        See Also
        --------
        scipy.stats.rankdata

        Examples
        --------
        >>> import numpy as np
        >>> from biom import Table
        >>> data = np.array([[ 99,  12,   8], [  0,  42,   7],
        ...                  [112,  42,   6], [  5,  75,   5]])
        >>> t = Table(data, sample_ids=['s1', 's2', 's3'],
        ...           observation_ids=['o1', 'o2', 'o3', 'o4'])

        Convert observation counts to their ranked abundance, from smallest
        to largest.

        >>> print(t.rankdata())  # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	s1	s2	s3
        o1	2.0	1.0	4.0
        o2	0.0	2.5	3.0
        o3	3.0	2.5	2.0
        o4	1.0	4.0	1.0

        """
        def f(val, id_, _):
            return scipy.stats.rankdata(val, method=method)
        return self.transform(f, axis=axis, inplace=inplace)

    def norm(self, axis='sample', inplace=True):
        """Normalize in place sample values by an observation, or vice versa.

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
            The axis to use for normalization.
        inplace : bool, optional
            Defaults to ``True``. If ``True``, performs the normalization in
            place. Otherwise, returns a new table with the normalization
            applied.

        Returns
        -------
        biom.Table
            The normalized table

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x2 table:

        >>> data = np.asarray([[2, 0], [6, 1]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2'])

        Get a version of the table normalized on the 'sample' axis, leaving the
        original table untouched:

        >>> new_table = table.norm(inplace=False)
        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2
        O1  2.0 0.0
        O2  6.0 1.0
        >>> print(new_table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2
        O1  0.25    0.0
        O2  0.75    1.0

        Get a version of the table normalized on the 'observation' axis,
        again leaving the original table untouched:

        >>> new_table = table.norm(axis='observation', inplace=False)
        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2
        O1  2.0 0.0
        O2  6.0 1.0
        >>> print(new_table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2
        O1  1.0 0.0
        O2  0.857142857143  0.142857142857

        Do the same normalization on 'observation', this time in-place:

        >>> table.norm(axis='observation')
        2 x 2 <class 'biom.table.Table'> with 3 nonzero entries (75% dense)
        >>> print(table) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID S1  S2
        O1  1.0 0.0
        O2  0.857142857143  0.142857142857
        """
        def f(val, id_, _):
            return val / float(val.sum())

        return self.transform(f, axis=axis, inplace=inplace)

    def nonzero(self):
        """Yields locations of nonzero elements within the data matrix

        Returns
        -------
        generator
            Yields ``(observation_id, sample_id)`` for each nonzero element
        """
        csr = self._data.tocsr()
        samp_ids = self.ids()
        obs_ids = self.ids(axis='observation')

        indptr = csr.indptr
        indices = csr.indices
        for row_idx in range(indptr.size - 1):
            start = indptr[row_idx]
            end = indptr[row_idx+1]

            obs_id = obs_ids[row_idx]
            for col_idx in indices[start:end]:
                yield (obs_id, samp_ids[col_idx])

    def nonzero_counts(self, axis, binary=True):
        """Get nonzero summaries about an axis

        Parameters
        ----------
        axis : {'sample', 'observation', 'whole'}
            The axis on which to count nonzero entries
        binary : bool, optional
            Defaults to ``True``. If ``True``, return number of nonzero
            entries. If ``False``, sum the values of the entries.

        Returns
        -------
        numpy.array
            Counts in index order to the axis
        """
        if binary:
            dtype = 'int'

            def op(x):
                return x.nonzero()[0].size
        else:
            dtype = self.dtype

            def op(x):
                return x.sum()

        if axis in ('sample', 'observation'):
            # can use np.bincount for CSMat or ScipySparse
            result = zeros(len(self.ids(axis=axis)), dtype=dtype)
            for idx, vals in enumerate(self.iter_data(axis=axis)):
                result[idx] = op(vals)
        else:
            result = zeros(1, dtype=dtype)
            for vals in self.iter_data():
                result[0] += op(vals)

        return result

    def _union_id_order(self, a, b):
        """Determines merge order for id lists A and B"""
        all_ids = list(a[:])
        all_ids.extend(b[:])
        new_order = {}
        idx = 0
        for id_ in all_ids:
            if id_ not in new_order:
                new_order[id_] = idx
                idx += 1
        return new_order

    def _intersect_id_order(self, a, b):
        """Determines the merge order for id lists A and B"""
        all_b = set(b[:])
        new_order = {}
        idx = 0
        for id_ in a:
            if id_ in all_b:
                new_order[id_] = idx
                idx += 1
        return new_order

    def remove_empty(self, axis='whole', inplace=True):
        """Remove empty samples or observations from the table

        Parameters
        ----------
        axis : {'whole', 'sample', 'observation'}, optional
            The axis on which to operate.
        inplace : bool, optional
            If ``True`` vectors are removed in ``self``; if ``False`` the
            vectors are removed in a new table is returned.

        Raises
        ------
        UnknownAxisError
            If the axis is not recognized.

        Returns
        -------
        Table
            A table object with the zero'd rows, or columns removed as
            specified by the `axis` parameter.
        """
        if axis not in ['sample', 'observation', 'whole']:
            raise UnknownAxisError(axis)

        if inplace:
            table = self
        else:
            table = self.copy()

        if axis == 'whole':
            axes = ['sample', 'observation']
        else:
            axes = [axis]

        for ax in axes:
            table.filter(table.ids(axis=ax)[table.sum(axis=ax) > 0], axis=ax)

        return table

    def align_to(self, other, axis='detect'):
        """Align self to other over a requested axis

        Parameters
        ----------
        other : biom.Table
            The table to align too
        axis : str, optional, {sample, observation, both, detect}
            If 'sample' or 'observation', align to that axis. If 'both', align
            both axes. If 'detect', align what can be aligned.

        Raises
        ------
        DisjointIDError
            If the requested axis can't be aligned.
        UnknownAxisError
            If an unrecognized axis is specified.

        Examples
        --------
        Align one table to another, for instance a table of 16S data to a table
        of metagenomic data. In this example, we're aligning the samples of the
        two tables.

        >>> from biom import Table
        >>> import numpy as np
        >>> amplicon = Table(np.array([[0, 1, 2], [3, 4, 5]]),
        ...                  ['Ecoli', 'Staphylococcus'],
        ...                  ['S1', 'S2', 'S3'])
        >>> metag = Table(np.array([[6, 7, 8], [9, 10, 11]]),
        ...               ['geneA', 'geneB'],
        ...               ['S3', 'S2', 'S1'])
        >>> amplicon = amplicon.align_to(metag)
        >>> print(amplicon)  # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	S3	S2	S1
        Ecoli	2.0	1.0	0.0
        Staphylococcus	5.0	4.0	3.0
        """
        self_o = set(self.ids(axis='observation'))
        self_s = set(self.ids())
        other_o = set(other.ids(axis='observation'))
        other_s = set(other.ids())

        alignable_o = self_o == other_o
        alignable_s = self_s == other_s

        if axis == 'both' and not (alignable_o and alignable_s):
            raise DisjointIDError("Cannot align both axes")
        elif axis == 'sample' and not alignable_s:
            raise DisjointIDError("Cannot align samples")
        elif axis == 'observation' and not alignable_o:
            raise DisjointIDError("Cannot align observations")
        elif axis == 'detect' and not (alignable_o or alignable_s):
            raise DisjointIDError("Neither axis appears alignable")

        if axis == 'both':
            order = ['observation', 'sample']
        elif axis == 'detect':
            order = []
            if alignable_s:
                order.append('sample')
            if alignable_o:
                order.append('observation')
        elif axis == 'sample':
            order = ['sample']
        elif axis == 'observation':
            order = ['observation']
        else:
            raise UnknownAxisError("Unrecognized axis: %s" % axis)

        table = self
        for aln_axis in order:
            table = table.sort_order(other.ids(axis=aln_axis),
                                     axis=aln_axis)

        return table

    def concat(self, others, axis='sample'):
        """Concatenate tables if axis is disjoint

        Parameters
        ----------
        others : iterable of biom.Table, or a single biom.Table instance
            Tables to concatenate
        axis : {'sample', 'observation'}, optional
            The axis to concatenate on. i.e., if axis is 'sample', then tables
            will be joined such that the set of sample IDs in the resulting
            table will be the union of sample IDs across all tables in others.

        Raises
        ------
        DisjointIDError
            If IDs over the axis are not disjoint.

        Notes
        -----
        The type of the table is inherited from self.

        Examples
        --------
        Concatenate three tables in which the sample IDs are disjoint. Note
        the observation IDs in this example are not disjoint (although they
        can be):

        >>> from biom import Table
        >>> import numpy as np
        >>> a = Table(np.array([[0, 1, 2], [3, 4, 5]]), ['O1', 'O2'],
        ...                     ['S1', 'S2', 'S3'],
        ...                     [{'taxonomy': 'foo'}, {'taxonomy': 'bar'}])
        >>> b = Table(np.array([[6, 7, 8], [9, 10, 11]]), ['O3', 'O4'],
        ...                     ['S4', 'S5', 'S6'],
        ...                     [{'taxonomy': 'baz'}, {'taxonomy': 'foobar'}])
        >>> c = Table(np.array([[12, 13, 14], [15, 16, 17]]), ['O1', 'O5'],
        ...                     ['S7', 'S8', 'S9'],
        ...                     [{'taxonomy': 'foo'}, {'taxonomy': 'biz'}])
        >>> d = a.concat([b, c])
        >>> print(d)  # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	S1	S2	S3	S4	S5	S6	S7	S8	S9
        O1	0.0	1.0	2.0	0.0	0.0	0.0	12.0	13.0	14.0
        O2	3.0	4.0	5.0	0.0	0.0	0.0	0.0	0.0	0.0
        O3	0.0	0.0	0.0	6.0	7.0	8.0	0.0	0.0	0.0
        O4	0.0	0.0	0.0	9.0	10.0	11.0	0.0	0.0	0.0
        O5	0.0	0.0	0.0	0.0	0.0	0.0	15.0	16.0	17.0

        """
        if isinstance(others, self.__class__):
            others = [others, ]

        # we grow along the opposite axis
        invaxis = self._invert_axis(axis)
        if axis == 'sample':
            dim_getter = itemgetter(1)
            stack = hstack
            invstack = vstack
        else:
            dim_getter = itemgetter(0)
            stack = vstack
            invstack = hstack

        axis_ids = set()
        invaxis_ids = set()
        invaxis_metadata = {}

        all_tables = others[:]
        all_tables.insert(0, self)

        # verify disjoint, and fetch all ids from all tables
        for table in all_tables:
            table_axis_ids = table.ids(axis=axis)
            table_invaxis_order = table.ids(axis=invaxis)
            table_invaxis = set(table_invaxis_order)

            # test we have disjoint IDs
            if not axis_ids.isdisjoint(table_axis_ids):
                raise DisjointIDError("IDs are not disjoint")
            axis_ids.update(table_axis_ids)

            if set(table_invaxis) - invaxis_ids:
                for i in (set(table_invaxis) - invaxis_ids):
                    invaxis_metadata[i] = table.metadata(i, axis=invaxis)

                # add to our perspective all inv axis IDs
                invaxis_ids.update(table_invaxis)

        invaxis_order = sorted(invaxis_ids)

        # determine what inv axis IDs do not exist per table, and pad and sort
        # as necessary
        padded_tables = []
        for table in all_tables:
            missing_ids = list(invaxis_ids - set(table.ids(axis=invaxis)))

            if missing_ids:
                # determine new shape
                n_invaxis = len(missing_ids)
                n_axis = len(table.ids(axis=axis))
                if axis == 'sample':
                    shape = (n_invaxis, n_axis)
                else:
                    shape = (n_axis, n_invaxis)

                # create the padded matrix
                zerod = csr_matrix(shape)
                tmp_mat = invstack([table.matrix_data, zerod])

                # resolve invert axis ids and metadata
                tmp_inv_ids = list(table.ids(axis=invaxis))
                tmp_inv_ids.extend(missing_ids)
                tmp_inv_md = table.metadata(axis=invaxis)
                if tmp_inv_md is None:
                    tmp_inv_md = [None] * len(table.ids(axis=invaxis))
                else:
                    tmp_inv_md = list(tmp_inv_md)
                tmp_inv_md.extend([invaxis_metadata[i] for i in missing_ids])

                # resolve axis ids and metadata
                tmp_ids = list(table.ids(axis=axis))
                tmp_md = table.metadata(axis=axis)

                # resolve construction based off axis. This really should be
                # pushed to a classmethod.
                if axis == 'sample':
                    tmp_table = self.__class__(tmp_mat, tmp_inv_ids, tmp_ids,
                                               tmp_inv_md, tmp_md)
                else:
                    tmp_table = self.__class__(tmp_mat, tmp_ids, tmp_inv_ids,
                                               tmp_md, tmp_inv_md)
            else:
                tmp_table = table

            # sort the table if necessary
            if (tmp_table.ids(axis=invaxis) == invaxis_order).all():
                padded_tables.append(tmp_table)
            else:
                padded_tables.append(tmp_table.sort_order(invaxis_order,
                                                          axis=invaxis))

        # actually concatenate the matrices, IDs and metadata
        concat_mat = stack([t.matrix_data for t in padded_tables])
        concat_ids = np.concatenate([t.ids(axis=axis) for t in padded_tables])
        concat_md = []
        for table in padded_tables:
            metadata = table.metadata(axis=axis)
            if metadata is None:
                metadata = [None] * dim_getter(table.shape)
            concat_md.extend(metadata)

        # inverse axis sourced from whatever is in the first table
        inv_md = padded_tables[0].metadata(axis=invaxis)
        if axis == 'sample':
            concat = self.__class__(concat_mat, invaxis_order, concat_ids,
                                    inv_md, concat_md, type=self.type)
        else:
            concat = self.__class__(concat_mat, concat_ids, invaxis_order,
                                    concat_md, inv_md, type=self.type)

        return concat

    def _fast_merge(self, others):
        """For simple merge operations it is faster to aggregate using pandas

        Parameters
        ----------
        others : Table, or Iterable of Table
            If a Table, then merge with that table. If an iterable, then merge
            all of the tables
        """
        tables = [self] + others

        # gather all identifiers across tables
        all_features = set(np.hstack([t.ids(axis='observation')
                                      for t in tables]))
        all_samples = set(np.hstack([t.ids() for t in tables]))

        # produce a new stable order
        feature_order = sorted(all_features)
        sample_order = sorted(all_samples)

        # generate unique integer ids for the identifiers, and let's order
        # it to be polite
        feature_map = {i: idx for idx, i in enumerate(feature_order)}
        sample_map = {i: idx for idx, i in enumerate(sample_order)}

        ntuples = sum([t.nnz for t in tables])

        # we're going to aggregate in COO. per scipy, it is efficient for
        # construction of large matrices. importantly, it allows for
        # duplicates which in this case correspond to multiple values for
        # the same sample/feature across tables. the duplicates are summed
        # implicitly on conversion to csr/csc.
        rows = np.empty(ntuples, dtype=np.int32)
        cols = np.empty(ntuples, dtype=np.int32)
        data = np.empty(ntuples, dtype=self.matrix_data.dtype)

        offset = 0
        for table in tables:
            t_nnz = table.nnz

            coo = table.matrix_data.tocoo()

            # we need to map the index positions in the current table to the
            # index positions in the full matrix
            row_map = np.array([feature_map[i]
                                for i in table.ids(axis='observation')],
                               dtype=np.int32)
            col_map = np.array([sample_map[i]
                                for i in table.ids()],
                               dtype=np.int32)
            coo.row = row_map[coo.row]
            coo.col = col_map[coo.col]

            # store our coo data
            rows[offset:offset + t_nnz] = coo.row
            cols[offset:offset + t_nnz] = coo.col
            data[offset:offset + t_nnz] = coo.data
            offset += t_nnz

        coo = coo_matrix((data, (rows, cols)),
                         shape=(len(feature_order), len(sample_order)))

        return self.__class__(coo.tocsr(), feature_order, sample_order)

    def merge(self, other, sample='union', observation='union',
              sample_metadata_f=prefer_self,
              observation_metadata_f=prefer_self):
        """Merge two tables together

        The axes, samples and observations, can be controlled independently.
        Both can work on either "union" or "intersection".

        `sample_metadata_f` and `observation_metadata_f` define how to
        merge metadata between tables. The default is to just keep the metadata
        associated to self if self has metadata otherwise take metadata from
        other. These functions are given both metadata dicts and must return
        a single metadata dict

        Parameters
        ----------
        other : biom.Table or Iterable of Table
            The other table to merge with this one. If an iterable, the tables
            are expected to not have metadata.
        sample : 'union', 'intersection', optional
            How the sample axis is handled
        observation : 'union', 'intersection', optional
            How the observation axis is handled
        sample_metadata_f : function, optional
            Defaults to ``biom.util.prefer_self``. Defines how to handle sample
            metadata during merge.
        obesrvation_metadata_f : function, optional
            Defaults to ``biom.util.prefer_self``. Defines how to handle
            observation metdata during merge.

        Returns
        -------
        biom.Table
            The merged table

        Notes
        -----
        - If ``sample_metadata_f`` and ``observation_metadata_f`` are None,
            then a fast merge is applied.
        - There is an implicit type conversion to ``float``.
        - The return type is always that of ``self``

        Examples
        --------

        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x2 table and a 3x2 table:

        >>> d_a = np.asarray([[2, 0], [6, 1]])
        >>> t_a = Table(d_a, ['O1', 'O2'], ['S1', 'S2'])
        >>> d_b = np.asarray([[4, 5], [0, 3], [10, 10]])
        >>> t_b = Table(d_b, ['O1', 'O2', 'O3'], ['S1', 'S2'])

        Merging the table results in the overlapping samples/observations (see
        `O1` and `S2`) to be summed and the non-overlapping ones to be added to
        the resulting table (see `S3`).

        >>> merged_table = t_a.merge(t_b)
        >>> print(merged_table)  # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	S1	S2
        O1	6.0	5.0
        O2	6.0	4.0
        O3	10.0	10.0

        """
        s_md = self.metadata()
        o_md = self.metadata(axis='observation')
        no_md = (s_md is None) and (o_md is None)
        ignore_md = (sample_metadata_f is None) and \
            (observation_metadata_f is None)

        if no_md or ignore_md:
            if sample == 'union' and observation == 'union':
                if isinstance(other, (list, set, tuple)):
                    return self._fast_merge(other)
                else:
                    return self._fast_merge([other, ])

        # determine the sample order in the resulting table
        if sample == 'union':
            new_samp_order = self._union_id_order(self.ids(), other.ids())
        elif sample == 'intersection':
            new_samp_order = self._intersect_id_order(self.ids(), other.ids())
        else:
            raise TableException("Unknown sample merge type: %s" % sample)

        # determine the observation order in the resulting table
        if observation == 'union':
            new_obs_order = self._union_id_order(
                self.ids(axis='observation'), other.ids(axis='observation'))
        elif observation == 'intersection':
            new_obs_order = self._intersect_id_order(
                self.ids(axis='observation'), other.ids(axis='observation'))
        else:
            raise TableException(
                "Unknown observation merge type: %s" %
                observation)

        # convert these to lists, no need to be dictionaries and reduces
        # calls to items() and allows for pre-caluculating insert order
        new_samp_order = sorted(new_samp_order.items(), key=itemgetter(1))
        new_obs_order = sorted(new_obs_order.items(), key=itemgetter(1))

        # if we don't have any samples, complain loudly. This is likely from
        # performing an intersection without overlapping ids
        if not new_samp_order:
            raise TableException("No samples in resulting table!")
        if not new_obs_order:
            raise TableException("No observations in resulting table!")

        # helper index lookups
        other_obs_idx = other._obs_index
        self_obs_idx = self._obs_index
        other_samp_idx = other._sample_index
        self_samp_idx = self._sample_index

        # pre-calculate sample order from each table. We only need to do this
        # once which dramatically reduces the number of dict lookups necessary
        # within the inner loop
        other_samp_order = []
        self_samp_order = []
        for samp_id, nsi in new_samp_order:  # nsi -> new_sample_index
            other_samp_order.append((nsi, other_samp_idx.get(samp_id, None)))
            self_samp_order.append((nsi, self_samp_idx.get(samp_id, None)))

        # pre-allocate the a list for placing the resulting vectors as the
        # placement id is not ordered
        vals = [None for i in range(len(new_obs_order))]

        # POSSIBLE DECOMPOSITION
        # resulting sample ids and sample metadata
        sample_ids = []
        sample_md = []
        self_sample_md = self.metadata()
        other_sample_md = other.metadata()
        for id_, idx in new_samp_order:
            sample_ids.append(id_)

            # if we have sample metadata, grab it
            if self_sample_md is None or not self.exists(id_):
                self_md = None
            else:
                self_md = self_sample_md[self_samp_idx[id_]]

            # if we have sample metadata, grab it
            if other_sample_md is None or not other.exists(id_):
                other_md = None
            else:
                other_md = other_sample_md[other_samp_idx[id_]]

            sample_md.append(sample_metadata_f(self_md, other_md))

        # POSSIBLE DECOMPOSITION
        # resulting observation ids and sample metadata
        obs_ids = []
        obs_md = []
        self_obs_md = self.metadata(axis='observation')
        other_obs_md = other.metadata(axis='observation')
        for id_, idx in new_obs_order:
            obs_ids.append(id_)

            # if we have observation metadata, grab it
            if self_obs_md is None or not self.exists(id_, axis="observation"):
                self_md = None
            else:
                self_md = self_obs_md[self_obs_idx[id_]]

            # if we have observation metadata, grab it
            if other_obs_md is None or \
                    not other.exists(id_, axis="observation"):
                other_md = None
            else:
                other_md = other_obs_md[other_obs_idx[id_]]

            obs_md.append(observation_metadata_f(self_md, other_md))

        # length used for construction of new vectors
        vec_length = len(new_samp_order)

        # walk over observations in our new order
        for obs_id, new_obs_idx in new_obs_order:
            # create new vector for matrix values
            new_vec = zeros(vec_length, dtype='float')

            # This method allows for the creation of a matrix of self type.
            # See note above
            # new_vec = data_f()

            # see if the observation exists in other, if so, pull it out.
            # if not, set to the placeholder missing
            if other.exists(obs_id, axis="observation"):
                other_vec = other.data(obs_id, 'observation')
            else:
                other_vec = None

            # see if the observation exists in self, if so, pull it out.
            # if not, set to the placeholder missing
            if self.exists(obs_id, axis="observation"):
                self_vec = self.data(obs_id, 'observation')
            else:
                self_vec = None

            # short circuit. If other doesn't have any values, then we can just
            # take all values from self
            if other_vec is None:
                for (n_idx, s_idx) in self_samp_order:
                    if s_idx is not None:
                        new_vec[n_idx] = self_vec[s_idx]

            # short circuit. If self doesn't have any values, then we can just
            # take all values from other
            elif self_vec is None:
                for (n_idx, o_idx) in other_samp_order:
                    if o_idx is not None:
                        new_vec[n_idx] = other_vec[o_idx]

            else:
                # NOTE: DM 7.5.12, no observed improvement at the profile level
                # was made on this inner loop by using self_samp_order and
                # other_samp_order lists.

                # walk over samples in our new order
                for samp_id, new_samp_idx in new_samp_order:
                    # pull out each individual sample value. This is expensive,
                    # but the vectors are in a different alignment. It is
                    # possible that this could be improved with numpy take but
                    # needs to handle missing values appropriately
                    if samp_id not in self_samp_idx:
                        self_vec_value = 0
                    else:
                        self_vec_value = self_vec[self_samp_idx[samp_id]]

                    if samp_id not in other_samp_idx:
                        other_vec_value = 0
                    else:
                        other_vec_value = other_vec[other_samp_idx[samp_id]]

                    new_vec[new_samp_idx] = self_vec_value + other_vec_value

            # convert our new vector to self type as to make sure we don't
            # accidently force a dense representation in memory
            vals[new_obs_idx] = self._conv_to_self_type(new_vec)

        return self.__class__(self._conv_to_self_type(vals), obs_ids[:],
                              sample_ids[:], obs_md, sample_md)

    @classmethod
    def from_hdf5(cls, h5grp, ids=None, axis='sample', parse_fs=None,
                  subset_with_metadata=True):
        """Parse an HDF5 formatted BIOM table

        If ids is provided, only the samples/observations listed in ids
        (depending on the value of axis) will be loaded

        The expected structure of this group is below. A few basic definitions,
        N is the number of observations and M is the number of samples. Data
        are stored in both compressed sparse row (for observation oriented
        operations) and compressed sparse column (for sample oriented
        operations).

        Notes
        -----
        The expected HDF5 group structure is below. An example of an HDF5 file
        in DDL can be found here [1]_.

        - ./id                                                  : str, an arbitrary ID  # noqa
        - ./type                                                : str, the table type (e.g, OTU table)  # noqa
        - ./format-url                                          : str, a URL that describes the format  # noqa
        - ./format-version                                      : two element tuple of int32, major and minor  # noqa
        - ./generated-by                                        : str, what generated this file  # noqa
        - ./creation-date                                       : str, ISO format  # noqa
        - ./shape                                               : two element tuple of int32, N by M  # noqa
        - ./nnz                                                 : int32 or int64, number of non zero elems  # noqa
        - ./observation                                         : Group  # noqa
        - ./observation/ids                                     : (N,) dataset of str or vlen str  # noqa
        - ./observation/matrix                                  : Group  # noqa
        - ./observation/matrix/data                             : (nnz,) dataset of float64  # noqa
        - ./observation/matrix/indices                          : (nnz,) dataset of int32  # noqa
        - ./observation/matrix/indptr                           : (M+1,) dataset of int32  # noqa
        - ./observation/metadata                                : Group  # noqa
        - [./observation/metadata/foo]                          : Optional, (N,) dataset of any valid HDF5 type in index order with IDs.  # noqa
        - ./observation/group-metadata                          : Group  # noqa
        - [./observation/group-metadata/foo]                    : Optional, (?,) dataset of group metadata that relates IDs  # noqa
        - [./observation/group-metadata/foo.attrs['data_type']] : attribute of the foo dataset that describes contained type (e.g., newick)  # noqa
        - ./sample                                              : Group  # noqa
        - ./sample/ids                                          : (M,) dataset of str or vlen str  # noqa
        - ./sample/matrix                                       : Group  # noqa
        - ./sample/matrix/data                                  : (nnz,) dataset of float64  # noqa
        - ./sample/matrix/indices                               : (nnz,) dataset of int32  # noqa
        - ./sample/matrix/indptr                                : (N+1,) dataset of int32  # noqa
        - ./sample/metadata                                     : Group  # noqa
        - [./sample/metadata/foo]                               : Optional, (M,) dataset of any valid HDF5 type in index order with IDs.  # noqa
        - ./sample/group-metadata                               : Group  # noqa
        - [./sample/group-metadata/foo]                         : Optional, (?,) dataset of group metadata that relates IDs  # noqa
        - [./sample/group-metadata/foo.attrs['data_type']]      : attribute of the foo dataset that describes contained type (e.g., newick)  # noqa

        The '?' character on the dataset size means that it can be of arbitrary
        length.

        The expected structure for each of the metadata datasets is a list of
        atomic type objects (int, float, str, ...), where the index order of
        the list corresponds to the index order of the relevant axis IDs.
        Special metadata fields have been defined, and they are stored in a
        specific way. Currently, the available special metadata fields are:

        - taxonomy: (N, ?) dataset of str or vlen str
        - KEGG_Pathways: (N, ?) dataset of str or vlen str
        - collapsed_ids: (N, ?) dataset of str or vlen str

        Parameters
        ----------
        h5grp : a h5py ``Group`` or an open h5py ``File``
            The object to load from
        ids : iterable
            The sample/observation ids of the samples/observations that we need
            to retrieve from the hdf5 biom table
        axis : 'sample', 'observation', optional
            The axis to subset on
        parse_fs : dict, optional
            Specify custom parsing functions for metadata fields. This dict is
            expected to be {'metadata_field': function}, where the function
            signature is (object) corresponding to a single row in the
            associated metadata dataset. The return from this function an
            object as well, and is the parsed representation of the metadata.
        subset_with_metadata : bool, optional
            When subsetting (i.e., `ids` is `not None`), whether to also parse
            the metadata. By default, the metadata are also subset. The reason
            for exposing this functionality is that, for large tables, there
            exists a very large overhead for this metadata manipulation.

        Returns
        -------
        biom.Table
            A BIOM ``Table`` object

        Raises
        ------
        ValueError
            If `ids` are not a subset of the samples or observations ids
            present in the hdf5 biom table
            If h5grp is not a HDF5 file or group

        References
        ----------
        .. [1] http://biom-format.org/documentation/format_versions/biom-2.1.\
html

        See Also
        --------
        Table.to_hdf5

        Examples
        --------
        >>> from biom.table import Table
        >>> from biom.util import biom_open
        >>> with biom_open('rich_sparse_otu_table_hdf5.biom') as f \
# doctest: +SKIP
        >>>     t = Table.from_hdf5(f) # doctest: +SKIP

        Parse a hdf5 biom table subsetting observations
        >>> from biom.util import biom_open # doctest: +SKIP
        >>> from biom.parse import parse_biom_table
        >>> with biom_open('rich_sparse_otu_table_hdf5.biom') as f \
# doctest: +SKIP
        >>>     t = Table.from_hdf5(f, ids=["GG_OTU_1"],
        ...                         axis='observation') # doctest: +SKIP
        """
        if not isinstance(h5grp, (h5py.Group, h5py.File)):
            raise ValueError("h5grp does not appear to be an HDF5 file or "
                             "group")

        if axis not in ['sample', 'observation']:
            raise UnknownAxisError(axis)

        if parse_fs is None:
            parse_fs = {}

        if not subset_with_metadata and ids is not None:
            ids = set(ids)

            raw_indices = h5grp['%s/matrix/indices' % axis]
            raw_indptr = h5grp['%s/matrix/indptr' % axis]
            raw_data = h5grp['%s/matrix/data' % axis]
            axis_ids = h5grp['%s/ids' % axis][:]

            to_keep = np.array([i for i, id_ in enumerate(axis_ids)
                                if id_ in ids], dtype=int)
            start_end = [(raw_indptr[i], raw_indptr[i+1]) for i in to_keep]
            indptr = np.empty(len(to_keep) + 1, dtype=np.int32)
            indptr[0] = 0
            indptr[1:] = np.array([e - s for s, e in start_end]).cumsum()
            data = np.concatenate([raw_data[s:e] for s, e in start_end])
            indices = np.concatenate([raw_indices[s:e] for s, e in start_end])

            if axis == 'sample':
                obs_ids = h5grp['observation/ids'][:]
                samp_ids = axis_ids[to_keep]
                shape = (len(obs_ids), len(to_keep))
                mat = csc_matrix((data, indices, indptr), shape=shape)
            else:
                samp_ids = h5grp['sample/ids'][:]
                obs_ids = axis_ids[to_keep]
                shape = (len(to_keep), len(samp_ids))
                mat = csr_matrix((data, indices, indptr), shape=shape)

            # use a fixed width dtype
            obs_ids_dtype = 'U%d' % max([len(v) for v in obs_ids])
            samp_ids_dtype = 'U%d' % max([len(v) for v in samp_ids])
            obs_ids = np.asarray(obs_ids, dtype=obs_ids_dtype)
            samp_ids = np.asarray(samp_ids, dtype=samp_ids_dtype)

            return Table(mat, obs_ids, samp_ids)

        id_ = h5grp.attrs['id']
        create_date = h5grp.attrs['creation-date']
        generated_by = h5grp.attrs['generated-by']

        if hasattr(datetime, "fromisoformat"):
            try:
                create_date = datetime.fromisoformat(create_date)
            except (TypeError, ValueError):
                pass

        shape = h5grp.attrs['shape']
        type_ = None if h5grp.attrs['type'] == '' else h5grp.attrs['type']

        if isinstance(id_, bytes):
            id_ = id_.decode('ascii')

        if isinstance(type_, bytes):
            type_ = type_.decode('ascii')

        def ensure_utf8(x):
            if isinstance(x, bytes):
                return x.decode('utf8')
            else:
                return

        def axis_load(grp):
            """Loads all the data of the given group"""
            # fetch all of the IDs
            ids = grp['ids'][:]

            if ids.size > 0:
                ids_dtype = 'U%d' % max([len(v) for v in ids])
                ids = np.asarray(ids, dtype=ids_dtype)

            parser = defaultdict(lambda: general_parser)
            parser['taxonomy'] = vlen_list_of_str_parser
            parser['Taxonomy'] = vlen_list_of_str_parser
            parser['KEGG_Pathways'] = vlen_list_of_str_parser
            parser['collapsed_ids'] = vlen_list_of_str_parser
            parser.update(parse_fs)

            # fetch ID specific metadata
            md = [{} for i in range(len(ids))]
            for category, dset in grp['metadata'].items():
                category = category.replace('@@SLASH@@', '/')
                parse_f = parser[category]
                data = dset[:]
                for md_dict, data_row in zip(md, data):
                    md_dict[category] = parse_f(data_row)

            # If there was no metadata on the axis, set it up as none
            md = md if any(md) else None

            # Fetch the group metadata
            grp_md = {cat: ensure_utf8(val[0])
                      for cat, val in grp['group-metadata'].items()}
            return ids, md, grp_md

        obs_ids, obs_md, obs_grp_md = axis_load(h5grp['observation'])
        samp_ids, samp_md, samp_grp_md = axis_load(h5grp['sample'])

        # load the data
        data_grp = h5grp[axis]['matrix']
        h5_data = data_grp["data"]
        h5_indices = data_grp["indices"]
        h5_indptr = data_grp["indptr"]

        # Check if we need to subset the biom table
        if ids is not None:
            def _get_ids(source_ids, desired_ids):
                """If desired_ids is not None, makes sure that it is a subset
                of source_ids and returns the desired_ids array-like and a
                boolean array indicating where the desired_ids can be found in
                source_ids"""
                if desired_ids is None:
                    ids = source_ids[:]
                    idx = np.ones(source_ids.shape, dtype=bool)
                else:
                    desired_ids = np.asarray(desired_ids)
                    # Get the index of the source ids to include
                    idx = np.in1d(source_ids, desired_ids)
                    # Retrieve only the ids that we are interested on
                    ids = source_ids[idx]
                    # Check that all desired ids have been found on source ids

                    if ids.shape != desired_ids.shape:
                        raise ValueError("The following ids could not be "
                                         "found in the biom table: %s" %
                                         (set(desired_ids) - set(ids)))
                return ids, idx

            # Get the observation and sample ids that we are interested in
            samp, obs = (ids, None) if axis == 'sample' else (None, ids)
            obs_ids, obs_idx = _get_ids(obs_ids, obs)
            samp_ids, samp_idx = _get_ids(samp_ids, samp)

            # Get the new matrix shape
            shape = (len(obs_ids), len(samp_ids))

            # Fetch the metadata that we are interested in
            def _subset_metadata(md, idx):
                """If md has data, returns the subset indicated by idx, a
                boolean array"""
                if md:
                    md = list(np.asarray(md)[np.where(idx)])
                return md

            obs_md = _subset_metadata(obs_md, obs_idx)
            samp_md = _subset_metadata(samp_md, samp_idx)

            # load the subset of the data
            idx = samp_idx if axis == 'sample' else obs_idx
            keep = np.where(idx)[0]
            indptr_indices = sorted(
                (h5_indptr[i], h5_indptr[i+1]) for i in keep
            )
            # Create the new indptr
            indptr_subset = np.array([end - start
                                      for start, end in indptr_indices])
            indptr = np.empty(len(keep) + 1, dtype=np.int32)
            indptr[0] = 0
            indptr[1:] = indptr_subset.cumsum()

            data = np.hstack([h5_data[start:end]
                              for start, end in indptr_indices])
            indices = np.hstack([h5_indices[start:end]
                                 for start, end in indptr_indices])
        else:
            # no subset need, just pass all data to scipy
            data = h5_data
            indices = h5_indices
            indptr = h5_indptr

        cs = (data, indices, indptr)

        if axis == 'sample':
            matrix = csc_matrix(cs, shape=shape)
        else:
            matrix = csr_matrix(cs, shape=shape)

        t = Table(matrix, obs_ids, samp_ids, obs_md or None,
                  samp_md or None, type=type_, create_date=create_date,
                  generated_by=generated_by, table_id=id_,
                  observation_group_metadata=obs_grp_md,
                  sample_group_metadata=samp_grp_md)

        if ids is not None:
            # filter out any empty samples or observations which may exist due
            # to subsetting
            def any_value(vals, id_, md):
                return np.any(vals)

            axis = 'observation' if axis == 'sample' else 'sample'
            t.filter(any_value, axis=axis)

        return t

    def to_dataframe(self, dense=False):
        """Convert matrix data to a Pandas SparseDataFrame or DataFrame

        Parameters
        ----------
        dense : bool, optional
            If True, return pd.DataFrame instead of pd.SparseDataFrame.

        Returns
        -------
        pd.DataFrame or pd.SparseDataFrame
            A DataFrame indexed on the observation IDs, with the column
            names as the sample IDs.

        Notes
        -----
        Metadata are not included.

        Examples
        --------
        >>> from biom import example_table
        >>> df = example_table.to_dataframe()
        >>> df
             S1   S2   S3
        O1  0.0  1.0  2.0
        O2  3.0  4.0  5.0
        """
        index = self.ids(axis='observation')
        columns = self.ids()

        if dense:
            mat = self.matrix_data.toarray()
            constructor = pd.DataFrame
        else:
            mat = self.matrix_data.copy()
            constructor = partial(pd.DataFrame.sparse.from_spmatrix)

        return constructor(mat, index=index, columns=columns)

    def to_anndata(self, dense=False, dtype="float32", transpose=True):
        """Convert Table to AnnData format

        Parameters
        ----------
        dense : bool, optional
            If True, set adata.X as np.ndarray instead of sparse matrix.
        dtype: str, optional
            dtype used for storage in anndata object.
        tranpose: bool, optional
            If True, transpose the anndata so that observations are columns

        Returns
        -------
        anndata.AnnData
            AnnData with matrix data and associated observation and
            sample metadata.

        Notes
        -----
        Nested metadata are not included.

        Examples
        --------
        >>> from biom import example_table
        >>> adata = example_table.to_anndata()
        >>> adata
        AnnData object with n_obs × n_vars = 3 × 2
            obs: 'environment'
            var: 'taxonomy_0', 'taxonomy_1'
        """
        try:
            import anndata
        except ImportError:
            raise ImportError(
                "Please install anndata package -- `pip install anndata`"
            )
        mat = self.matrix_data

        if dense:
            mat = mat.toarray()

        var = self.metadata_to_dataframe("sample")
        obs = self.metadata_to_dataframe("observation")

        adata = anndata.AnnData(mat, obs=obs, var=var, dtype=dtype)
        # Convention for scRNA-seq analysis in Python
        adata = adata.transpose()

        return adata

    def metadata_to_dataframe(self, axis):
        """Convert axis metadata to a Pandas DataFrame

        Parameters
        ----------
        axis : {'sample', 'observation'}
            The axis to operate on.

        Returns
        -------
        pd.DataFrame
            A DataFrame indexed by the ids of the desired axis, columns by the
            metadata keys over that axis.

        Raises
        ------
        UnknownAxisError
            If the requested axis isn't recognized
        KeyError
            IF the requested axis does not have metadata
        TypeError
            If a metadata column is a list or tuple, but is jagged over the
            axis.

        Notes
        -----
        Nested metadata (e.g., KEGG_Pathways) is not supported.

        Metadata which are lists or tuples (e.g., taxonomy) are expanded such
        that each index position is a unique column. For instance, the key
        taxonomy will become "taxonomy_0", "taxonomy_1", etc where "taxonomy_0"
        corresponds to the 0th index position of the taxonomy.

        Examples
        --------
        >>> from biom import example_table
        >>> example_table.metadata_to_dataframe('observation')
           taxonomy_0     taxonomy_1
        O1   Bacteria     Firmicutes
        O2   Bacteria  Bacteroidetes
        """
        md = self.metadata(axis=axis)
        if md is None:
            raise KeyError("%s does not have metadata" % axis)

        mcols = []
        for test in md:
            columns = []
            expand = {}
            for key, value in test.items():
                if isinstance(value, (tuple, list)):
                    expand[key] = True
                    for idx in range(len(value)):
                        columns.append("%s_%d" % (key, idx))
                else:
                    expand[key] = False
                    columns.append(key)
            if len(columns) > len(mcols):
                mcols = columns

        rows = []
        for m in md:
            row = []
            for key, value in m.items():
                if expand[key]:
                    for v in value:
                        row.append(v)
                else:
                    row.append(value)
            rows.append(row)

        return pd.DataFrame(rows, index=self.ids(axis=axis), columns=mcols)

    def to_hdf5(self, h5grp, generated_by, compress=True, format_fs=None,
                creation_date=None):
        """Store CSC and CSR in place

        The resulting structure of this group is below. A few basic
        definitions, N is the number of observations and M is the number of
        samples. Data are stored in both compressed sparse row [1]_ (CSR, for
        observation oriented operations) and compressed sparse column [2]_
        (CSC, for sample oriented operations).

        Notes
        -----
        This method does not return anything and operates in place on h5grp.

        The expected HDF5 group structure is below. An example of an HDF5 file
        in DDL can be found here [3]_.

        - ./id                                                  : str, an \
arbitrary ID
        - ./type                                                : str, the \
table type (e.g, OTU table)
        - ./format-url                                          : str, a URL \
that describes the format
        - ./format-version                                      : two element \
tuple of int32, major and minor
        - ./generated-by                                        : str, what \
generated this file
        - ./creation-date                                       : str, ISO \
format
        - ./shape                                               : two element \
tuple of int32, N by M
        - ./nnz                                                 : int32 or \
int64, number of non zero elems
        - ./observation                                         : Group
        - ./observation/ids                                     : (N,) dataset\
 of str or vlen str
        - ./observation/matrix                                  : Group
        - ./observation/matrix/data                             : (nnz,) \
dataset of float64
        - ./observation/matrix/indices                          : (nnz,) \
dataset of int32
        - ./observation/matrix/indptr                           : (M+1,) \
dataset of int32
        - ./observation/metadata                                : Group
        - [./observation/metadata/foo]                          : Optional, \
(N,) dataset of any valid HDF5 type in index order with IDs.
        - ./observation/group-metadata                          : Group
        - [./observation/group-metadata/foo]                    : Optional, \
(?,) dataset of group metadata that relates IDs
        - [./observation/group-metadata/foo.attrs['data_type']] : attribute of\
 the foo dataset that describes contained type (e.g., newick)
        - ./sample                                              : Group
        - ./sample/ids                                          : (M,) dataset\
 of str or vlen str
        - ./sample/matrix                                       : Group
        - ./sample/matrix/data                                  : (nnz,) \
dataset of float64
        - ./sample/matrix/indices                               : (nnz,) \
dataset of int32
        - ./sample/matrix/indptr                                : (N+1,) \
dataset of int32
        - ./sample/metadata                                     : Group
        - [./sample/metadata/foo]                               : Optional, \
(M,) dataset of any valid HDF5 type in index order with IDs.
        - ./sample/group-metadata                               : Group
        - [./sample/group-metadata/foo]                         : Optional, \
(?,) dataset of group metadata that relates IDs
        - [./sample/group-metadata/foo.attrs['data_type']]      : attribute of\
 the foo dataset that describes contained type (e.g., newick)

        The '?' character on the dataset size means that it can be of arbitrary
        length.

        The expected structure for each of the metadata datasets is a list of
        atomic type objects (int, float, str, ...), where the index order of
        the list corresponds to the index order of the relevant axis IDs.
        Special metadata fields have been defined, and they are stored in a
        specific way. Currently, the available special metadata fields are:

        - taxonomy: (N, ?) dataset of str or vlen str
        - KEGG_Pathways: (N, ?) dataset of str or vlen str
        - collapsed_ids: (N, ?) dataset of str or vlen str

        Parameters
        ----------
        h5grp : `h5py.Group` or `h5py.File`
            The HDF5 entity in which to write the BIOM formatted data.
        generated_by : str
            A description of what generated the table
        compress : bool, optional
            Defaults to ``True`` means fields will be compressed with gzip,
            ``False`` means no compression
        format_fs : dict, optional
            Specify custom formatting functions for metadata fields. This dict
            is expected to be {'metadata_field': function}, where the function
            signature is (h5py.Group, str, dict, bool) corresponding to the
            specific HDF5 group the metadata dataset will be associated with,
            the category being operated on, the metadata for the entire axis
            being operated on, and whether to enable compression on the
            dataset.  Anything returned by this function is ignored.
        creation_date : datetime, optional
            If provided, use this specific datetime on write as the creation
            timestamp

        See Also
        --------
        Table.from_hdf5

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/sci\
py.sparse.csr_matrix.html
        .. [2] http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/sci\
py.sparse.csc_matrix.html
        .. [3] http://biom-format.org/documentation/format_versions/biom-2.1.\
html

        Examples
        --------
        >>> from biom.util import biom_open  # doctest: +SKIP
        >>> from biom.table import Table
        >>> from numpy import array
        >>> t = Table(array([[1, 2], [3, 4]]), ['a', 'b'], ['x', 'y'])
        >>> with biom_open('foo.biom', 'w') as f:  # doctest: +SKIP
        ...     t.to_hdf5(f, "example")

        """
        if format_fs is None:
            format_fs = {}

        # cache nnz to avoid multiple calls to eliminate_zeros() as it is
        # expensive per profiling
        nnz = self.nnz
        h5grp.attrs['id'] = self.table_id if self.table_id else "No Table ID"
        h5grp.attrs['type'] = self.type if self.type else ""
        h5grp.attrs['format-url'] = "http://biom-format.org"
        h5grp.attrs['format-version'] = self.format_version
        h5grp.attrs['generated-by'] = generated_by
        if creation_date is None:
            h5grp.attrs['creation-date'] = datetime.now().isoformat()
        else:
            h5grp.attrs['creation-date'] = creation_date.isoformat()
        h5grp.attrs['shape'] = self.shape
        h5grp.attrs['nnz'] = nnz

        compression = None
        if compress is True:
            compression = 'gzip'

        formatter = defaultdict(lambda: general_formatter)
        formatter['taxonomy'] = vlen_list_of_str_formatter
        formatter['Taxonomy'] = vlen_list_of_str_formatter
        formatter['KEGG_Pathways'] = vlen_list_of_str_formatter
        formatter['collapsed_ids'] = vlen_list_of_str_formatter
        formatter.update(format_fs)

        for axis, order in zip(['observation', 'sample'], ['csr', 'csc']):
            grp = h5grp.create_group(axis)

            self._data = self._data.asformat(order)

            ids = self.ids(axis=axis)
            len_ids = len(ids)
            len_indptr = len(self._data.indptr)
            len_data = nnz

            md = self.metadata(axis=axis)

            # Create the group for the metadata
            grp.create_group('metadata')
            if md:
                exp = set(md[0])
                for other_id, other_md in zip(ids[1:], md[1:]):
                    if set(other_md) != exp:
                        raise ValueError("%s has inconsistent metadata "
                                         "categories with %s:\n"
                                         "%s: %s\n"
                                         "%s: %s" % (other_id, ids[0],
                                                     other_id, list(other_md),
                                                     ids[0], list(exp)))

                for category in list(md[0]):
                    # Create the dataset for the current category,
                    # putting values in id order
                    formatter[category](grp, category, md, compression)

            group_md = self.group_metadata(axis)

            # Create the group for the group metadata
            grp.create_group('group-metadata')

            if group_md:
                for key, value in group_md.items():
                    datatype, val = value
                    grp_dataset = grp.create_dataset(
                        'group-metadata/%s' % key,
                        shape=(1,), dtype=H5PY_VLEN_STR,
                        data=val, compression=compression)
                    grp_dataset.attrs['data_type'] = datatype

            grp.create_group('matrix')
            grp.create_dataset('matrix/data', shape=(len_data,),
                               dtype=np.float64,
                               data=self._data.data,
                               compression=compression)
            grp.create_dataset('matrix/indices', shape=(len_data,),
                               dtype=np.int32,
                               data=self._data.indices,
                               compression=compression)
            grp.create_dataset('matrix/indptr', shape=(len_indptr,),
                               dtype=np.int32,
                               data=self._data.indptr,
                               compression=compression)

            if len_ids > 0:
                # if we store IDs in the table as numpy arrays then this store
                # is cleaner, as is the parse
                grp.create_dataset('ids', shape=(len_ids,),
                                   dtype=H5PY_VLEN_STR,
                                   data=[i.encode('utf8') for i in ids],
                                   compression=compression)
            else:
                # Empty H5PY_VLEN_STR datasets are not supported.
                grp.create_dataset('ids', shape=(0, ), data=[],
                                   compression=compression)

    @classmethod
    def from_json(self, json_table, data_pump=None,
                  input_is_dense=False):
        """Parse a biom otu table type

        Parameters
        ----------
        json_table : dict
            A JSON object or dict that represents the BIOM table
        data_pump : tuple or None
            A secondary source of data
        input_is_dense : bool
            If `True`, the data contained will be interpretted as dense

        Returns
        -------
        Table

        Examples
        --------
        >>> from biom import Table
        >>> json_obj = {"id": "None",
        ...             "format": "Biological Observation Matrix 1.0.0",
        ...             "format_url": "http://biom-format.org",
        ...             "generated_by": "foo",
        ...             "type": "OTU table",
        ...             "date": "2014-06-03T14:24:40.884420",
        ...             "matrix_element_type": "float",
        ...             "shape": [5, 6],
        ...             "data": [[0,2,1.0],
        ...                      [1,0,5.0],
        ...                      [1,1,1.0],
        ...                      [1,3,2.0],
        ...                      [1,4,3.0],
        ...                      [1,5,1.0],
        ...                      [2,2,1.0],
        ...                      [2,3,4.0],
        ...                      [2,5,2.0],
        ...                      [3,0,2.0],
        ...                      [3,1,1.0],
        ...                      [3,2,1.0],
        ...                      [3,5,1.0],
        ...                      [4,1,1.0],
        ...                      [4,2,1.0]],
        ...             "rows": [{"id": "GG_OTU_1", "metadata": None},
        ...                      {"id": "GG_OTU_2", "metadata": None},
        ...                      {"id": "GG_OTU_3", "metadata": None},
        ...                      {"id": "GG_OTU_4", "metadata": None},
        ...                      {"id": "GG_OTU_5", "metadata": None}],
        ...             "columns": [{"id": "Sample1", "metadata": None},
        ...                         {"id": "Sample2", "metadata": None},
        ...                         {"id": "Sample3", "metadata": None},
        ...                         {"id": "Sample4", "metadata": None},
        ...                         {"id": "Sample5", "metadata": None},
        ...                         {"id": "Sample6", "metadata": None}]
        ...             }
        >>> t = Table.from_json(json_obj)

        """
        sample_ids = [col['id'] for col in json_table['columns']]
        sample_metadata = [col['metadata'] for col in json_table['columns']]
        obs_ids = [row['id'] for row in json_table['rows']]
        obs_metadata = [row['metadata'] for row in json_table['rows']]
        dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]
        if 'matrix_type' in json_table:
            if json_table['matrix_type'] == 'dense':
                input_is_dense = True
            else:
                input_is_dense = False
        type_ = json_table['type']

        if data_pump is None:
            data = json_table['data']
        else:
            data = data_pump
        create_date = None
        if hasattr(datetime, "fromisoformat"):
            try:
                create_date = datetime.fromisoformat(json_table['date'])
            except (TypeError, ValueError):
                pass
        table_obj = Table(data, obs_ids, sample_ids,
                          obs_metadata, sample_metadata,
                          shape=json_table['shape'],
                          dtype=dtype,
                          type=type_,
                          create_date=create_date,
                          generated_by=json_table['generated_by'],
                          input_is_dense=input_is_dense)
        return table_obj

    def to_json(self, generated_by, direct_io=None, creation_date=None):
        """Returns a JSON string representing the table in BIOM format.

        Parameters
        ----------
        generated_by : str
            a string describing the software used to build the table
        direct_io : file or file-like object, optional
            Defaults to ``None``. Must implementing a ``write`` function. If
            `direct_io` is not ``None``, the final output is written directly
            to `direct_io` during processing.
        creation_date : datetime, optional
            If provided, use this datetime as the creation date on write.

        Returns
        -------
        str
            A JSON-formatted string representing the biom table
        """
        if not isinstance(generated_by, str):
            raise TableException("Must specify a generated_by string")

        if creation_date is None:
            creation_date = datetime.now().isoformat()
        else:
            creation_date = creation_date.isoformat()

        # Fill in top-level metadata.
        if direct_io:
            direct_io.write('{')
            direct_io.write('"id": "%s",' % str(self.table_id))
            direct_io.write(
                '"format": "%s",' %
                get_biom_format_version_string((1, 0)))  # JSON table -> 1.0.0
            direct_io.write(
                '"format_url": "%s",' %
                get_biom_format_url_string())
            direct_io.write('"generated_by": "%s",' % generated_by)
            direct_io.write('"date": "%s",' % creation_date)
        else:
            id_ = '"id": "%s",' % str(self.table_id)
            format_ = '"format": "%s",' % get_biom_format_version_string(
                (1, 0))  # JSON table -> 1.0.0
            format_url = '"format_url": "%s",' % get_biom_format_url_string()
            generated_by = '"generated_by": "%s",' % generated_by
            date = '"date": "%s",' % creation_date

        # Determine if we have any data in the matrix, and what the shape of
        # the matrix is.
        try:
            num_rows, num_cols = self.shape
        except:  # noqa
            num_rows = num_cols = 0
        has_data = True if num_rows > 0 and num_cols > 0 else False

        # Default the matrix element type to test to be an integer in case we
        # don't have any data in the matrix to test.
        test_element = 0
        if has_data:
            test_element = self[0, 0]

        # Determine the type of elements the matrix is storing.
        if isinstance(test_element, int):
            matrix_element_type = "int"
        elif isinstance(test_element, float):
            matrix_element_type = "float"
        elif isinstance(test_element, str):
            matrix_element_type = "str"
        else:
            raise TableException("Unsupported matrix data type.")

        # Fill in details about the matrix.
        if direct_io:
            direct_io.write(
                '"matrix_element_type": "%s",' %
                matrix_element_type)
            direct_io.write('"shape": [%d, %d],' % (num_rows, num_cols))
        else:
            matrix_element_type = '"matrix_element_type": "%s",' % \
                matrix_element_type
            shape = '"shape": [%d, %d],' % (num_rows, num_cols)

        # Fill in the table type
        if self.type is None:
            type_ = '"type": null,'
        else:
            type_ = '"type": "%s",' % self.type

        if direct_io:
            direct_io.write(type_)

        # Fill in details about the rows in the table and fill in the matrix's
        # data. BIOM 2.0+ is now only sparse
        if direct_io:
            direct_io.write('"matrix_type": "sparse",')
            direct_io.write('"data": [')
        else:
            matrix_type = '"matrix_type": "sparse",'
            data = ['"data": [']

        max_row_idx = len(self.ids(axis='observation')) - 1
        max_col_idx = len(self.ids()) - 1
        rows = ['"rows": [']
        have_written = False
        for obs_index, obs in enumerate(self.iter(axis='observation')):
            # i'm crying on the inside
            if obs_index != max_row_idx:
                rows.append(
                    f'{{"id": {dumps(obs[1])}, "metadata": {dumps(obs[2])}}},'
                )
            else:
                rows.append(
                    f'{{"id": {dumps(obs[1])}, "metadata": {dumps(obs[2])}}}],'
                )

            # turns out its a pain to figure out when to place commas. the
            # simple work around, at the expense of a little memory
            # (bound by the number of samples) is to build of what will be
            # written, and then add in the commas where necessary.
            built_row = []
            for col_index, val in enumerate(obs[0]):
                if float(val) != 0.0:
                    built_row.append(
                        "[%d,%d,%f]" % (obs_index, col_index, val)
                    )
            if built_row:
                # if we have written a row already, its safe to add a comma
                if have_written:
                    if direct_io:
                        direct_io.write(',')
                    else:
                        data.append(',')
                if direct_io:
                    direct_io.write(','.join(built_row))
                else:
                    data.append(','.join(built_row))

                have_written = True

        # finalize the data block
        if direct_io:
            direct_io.write("],")
        else:
            data.append("],")

        # Fill in details about the columns in the table.
        columns = ['"columns": [']
        for samp_index, samp in enumerate(self.iter()):
            if samp_index != max_col_idx:
                columns.append('{{"id": {}, "metadata": {}}},'.format(
                    dumps(samp[1]), dumps(samp[2])))
            else:
                columns.append('{{"id": {}, "metadata": {}}}]'.format(
                    dumps(samp[1]), dumps(samp[2])))

        if rows[0] == '"rows": [' and len(rows) == 1:
            # empty table case
            rows = ['"rows": [],']
            columns = ['"columns": []']

        rows = ''.join(rows)
        columns = ''.join(columns)

        if direct_io:
            direct_io.write(rows)
            direct_io.write(columns)
            direct_io.write('}')
        else:
            return "{%s}" % ''.join([
                id_,
                format_,
                format_url,
                matrix_type,
                generated_by,
                date,
                type_,
                matrix_element_type,
                shape,
                ''.join(data),
                rows,
                columns,
            ])

    @staticmethod
    def from_adjacency(lines):
        """Parse an adjacency format into BIOM

        Parameters
        ----------
        lines : list, str, or file-like object
            The tab delimited data to parse

        Returns
        -------
        biom.Table
            A BIOM ``Table`` object

        Notes
        -----
        The input is expected to be of the form: observation, sample, value. A
        header is not required, but if present, it must be of the form:

        #OTU ID<tab>SampleID<tab>value

        Raises
        ------
        ValueError
            If the input is not an iterable or file-like object.
        ValueError
            If the data is incorrectly formatted.

        Examples
        --------
        Parse tab separated adjacency data into a table:

        >>> from biom.table import Table
        >>> from io import StringIO
        >>> data = 'a\\tb\\t1\\na\\tc\\t2\\nd\\tc\\t3'
        >>> data_fh = StringIO(data)
        >>> test_table = Table.from_adjacency(data_fh)
        """
        if not isinstance(lines, (list, tuple)):
            if hasattr(lines, 'readlines'):
                lines = lines.readlines()
            elif hasattr(lines, 'splitlines'):
                lines = lines.splitlines()
            else:
                raise ValueError("Not sure how to handle this input")

        def is_num(item):
            # from https://stackoverflow.com/a/23059703
            numeric = re.compile(r'(?=.)([+-]?([0-9]*)(\.([0-9]+))?)([eE][+-]?\d+)?')  # noqa
            match = numeric.match(item)
            start, stop = match.span()
            if (stop - start) == len(item):
                return True
            else:
                return False

        # sanity check and determine if we have a header or not
        lh = lines[0].strip().split('\t')
        if len(lh) != 3:
            raise ValueError("Does not appear to be an adjacency format")
        elif lh == ['#OTU ID', 'SampleID', 'value']:
            include_line_zero = False
        elif is_num(lh[2]):
            # allow anything for columns 1 and 2, but test that column 3 is
            # numeric
            include_line_zero = True
        else:
            raise ValueError("Does not appear to be an adjacency format")

        if not include_line_zero:
            lines = lines[1:]

        # extract the entities
        observations = []
        samples = []
        values = []
        for line in lines:
            parts = line.split('\t')
            assert len(parts) == 3
            observations.append(parts[0])
            samples.append(parts[1])
            values.append(float(parts[2]))

        # determine a stable order and index positioning for the identifiers
        obs_order = sorted(set(observations))
        samp_order = sorted(set(samples))
        obs_index = {o: i for i, o in enumerate(obs_order)}
        samp_index = {s: i for i, s in enumerate(samp_order)}

        # fill the matrix
        row = np.array([obs_index[obs] for obs in observations], dtype=int)
        col = np.array([samp_index[samp] for samp in samples], dtype=int)
        data = np.asarray(values)
        mat = coo_matrix((data, (row, col)))

        return Table(mat, obs_order, samp_order)

    @staticmethod
    def from_tsv(lines, obs_mapping, sample_mapping,
                 process_func, **kwargs):
        """Parse a tab separated (observation x sample) formatted BIOM table

        Parameters
        ----------
        lines : list, or file-like object
            The tab delimited data to parse
        obs_mapping : dict or None
            The corresponding observation metadata
        sample_mapping : dict or None
            The corresponding sample metadata
        process_func : function
            A function to transform the observation metadata

        Returns
        -------
        biom.Table
            A BIOM ``Table`` object

        Examples
        --------
        Parse tab separated data into a table:

        >>> from biom.table import Table
        >>> from io import StringIO
        >>> tsv = 'a\\tb\\tc\\n1\\t2\\t3\\n4\\t5\\t6'
        >>> tsv_fh = StringIO(tsv)
        >>> func = lambda x : x
        >>> test_table = Table.from_tsv(tsv_fh, None, None, func)
        """
        (sample_ids, obs_ids, data, t_md,
         t_md_name) = Table._extract_data_from_tsv(lines, **kwargs)

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

        return Table(data, obs_ids, sample_ids, obs_metadata, sample_metadata)

    @staticmethod
    def _extract_data_from_tsv(lines, delim='\t', dtype=float, md_parse=None):
        """Parse a classic table into (sample_ids, obs_ids, data, metadata,
        name)

        Parameters
        ----------
        lines: list or file-like object
            delimted data to parse
        delim: string
            delimeter in file lines
        dtype: type
            The expected type
        md_parse:  function or None
            funtion used to parse metdata

        Returns
        -------
        list
            sample_ids
        list
            observation_ids
        array
            data
        list
            metadata
        string
            column name if last column is non-numeric

        Notes
        ------
        This is intended to be close to how QIIME classic OTU tables are parsed
        with the exception of the additional md_name field

        This function is ported from QIIME (http://www.qiime.org), previously
        named parse_classic_otu_table. QIIME is a GPL project, but we obtained
        permission from the authors of this function to port it to the BIOM
        Format project (and keep it under BIOM's BSD license).

        .. shownumpydoc
        """
        def isfloat(value):
            # see https://stackoverflow.com/a/20929881
            try:
                float(value)
                return True
            except ValueError:
                return False

        if not isinstance(lines, list):
            try:
                hasattr(lines, 'seek')
            except AttributeError:
                raise RuntimeError(
                    "Input needs to support seek or be indexable")

        # find header, the first line that is not empty and does not start
        # with a #
        header = False
        list_index = 0
        data_start = 0
        for line in lines:
            if not line.strip():
                continue
            if not line.startswith('#'):
                # Covers the case where the first line is the header
                # and there is no indication of it (no comment character)
                if not header:
                    header = line.rstrip().split(delim)[1:]
                    data_start = list_index + 1
                else:
                    data_start = list_index
                break
            list_index += 1
            header = line.strip().split(delim)[1:]

        # If the first line is the header, then we need to get the data lines
        # line for the "last column" check
        if isinstance(lines, list):
            value_checks = lines[data_start:]
        else:
            lines.seek(0)
            for index in range(0, data_start):
                lines.readline()
            value_checks = [line for line in lines]

        # attempt to determine if the last column is non-numeric, ie, metadata
        last_values = [line.rsplit(delim, 1)[-1].strip()
                       for line in value_checks]
        last_column_is_numeric = all([isfloat(i) for i in last_values])

        # determine sample ids
        if last_column_is_numeric or data_start == 0:
            md_name = None
            metadata = None
            samp_ids = header[:]
        else:
            md_name = header[-1]
            metadata = []
            samp_ids = header[:-1]

        data = []
        obs_ids = []
        row_number = 0

        # Go back to the beginning if it is a file:
        if hasattr(lines, 'seek'):
            lines.seek(0)
            for index in range(0, data_start):
                line = lines.readline()
        else:
            lines = lines[data_start:]

        for lineno, line in enumerate(lines, data_start):
            if not line.strip():
                continue
            if line.startswith('#'):
                continue

            fields = line.split(delim)
            fields[-1] = fields[-1].strip()
            obs_ids.append(fields[0])

            if last_column_is_numeric:
                try:
                    values = list(map(dtype, fields[1:]))
                except ValueError:
                    badval, badidx = _identify_bad_value(dtype, fields[1:])
                    msg = "Invalid value on line %d, column %d, value %s"
                    raise TypeError(msg % (lineno, badidx+1, badval))
            else:
                try:
                    values = list(map(dtype, fields[1:-1]))
                except ValueError:
                    badval, badidx = _identify_bad_value(dtype, fields[1:])
                    msg = "Invalid value on line %d, column %d, value %s"
                    raise TypeError(msg % (lineno, badidx+1, badval))

                if md_parse is not None:
                    metadata.append(md_parse(fields[-1]))
                else:
                    metadata.append(fields[-1])
            for column_number in range(0, len(values)):
                if values[column_number] != dtype(0):
                    data.append([row_number, column_number,
                                 values[column_number]])
            row_number += 1
        return samp_ids, obs_ids, data, metadata, md_name

    def to_tsv(self, header_key=None, header_value=None,
               metadata_formatter=str,
               observation_column_name='#OTU ID',
               direct_io=None):
        """Return self as a string in tab delimited form

        Default ``str`` output for the ``Table`` is just row/col ids and table
        data without any metadata

        Parameters
        ----------
        header_key : str or ``None``, optional
            Defaults to ``None``
        header_value : str or ``None``, optional
            Defaults to ``None``
        metadata_formatter : function, optional
            Defaults to ``str``.  a function which takes a metadata entry and
            returns a formatted version that should be written to file
        observation_column_name : str, optional
            Defaults to "#OTU ID". The name of the first column in the output
            table, corresponding to the observation IDs.
        direct_io : file or file-like object, optional
            Defaults to ``None``. Must implement a ``write`` function. If
            `direct_io` is not ``None``, the final output is written directly
            to `direct_io` during processing.

        Returns
        -------
        str
            tab delimited representation of the Table

        Examples
        --------

        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table, with observation metadata and no sample
        metadata:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'],
        ...               [{'foo': 'bar'}, {'x': 'y'}], None)
        >>> print(table.to_tsv()) # doctest: +NORMALIZE_WHITESPACE
        # Constructed from biom file
        #OTU ID	S1	S2	S3
        O1	0.0	0.0	1.0
        O2	1.0	3.0	42.0
        >>> with open("result.tsv", "w") as f:
                table.to_tsv(direct_io=f)
        """
        return self.delimited_self('\t', header_key, header_value,
                                   metadata_formatter,
                                   observation_column_name,
                                   direct_io=direct_io)


def coo_arrays_to_sparse(data, dtype=np.float64, shape=None):
    """Map directly on to the coo_matrix constructor

    Parameters
    ----------
    data : tuple
        data must be (values, (rows, cols))
    dtype : type, optional
        Defaults to ``np.float64``
    shape : tuple or ``None``, optional
        Defaults to ``None``. If `shape` is ``None``, shape will be determined
        automatically from `data`.
    """
    if shape is None:
        values, (rows, cols) = data
        n_rows = max(rows) + 1
        n_cols = max(cols) + 1
    else:
        n_rows, n_cols = shape

    # coo_matrix allows zeros to be added as data, and this affects
    # nnz, items, and iteritems. Clean them out here, as this is
    # the only time these zeros can creep in.
    # Note: coo_matrix allows duplicate entries; the entries will
    # be summed when converted. Not really sure how we want to
    # handle this generally within BIOM- I'm okay with leaving it
    # as undefined behavior for now.
    matrix = coo_matrix(data, shape=(n_rows, n_cols), dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_list_to_sparse(data, dtype=float, shape=None):
    """Convert a list of lists into a scipy.sparse matrix.

    Parameters
    ----------
    data : iterable of iterables
        `data` should be in the format [[row, col, value], ...]
    dtype : type, optional
        defaults to ``float``
    shape : tuple or ``None``, optional
        Defaults to ``None``. If `shape` is ``None``, shape will be determined
        automatically from `data`.

    Returns
    -------
    scipy.csr_matrix
        The newly generated matrix
    """
    rows, cols, values = zip(*data)

    if shape is None:
        n_rows = max(rows) + 1
        n_cols = max(cols) + 1
    else:
        n_rows, n_cols = shape

    matrix = coo_matrix((values, (rows, cols)), shape=(n_rows, n_cols),
                        dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def nparray_to_sparse(data, dtype=float):
    """Convert a numpy array to a scipy.sparse matrix.

    Parameters
    ----------
    data : numpy.array
        The data to convert into a sparse matrix
    dtype : type, optional
        Defaults to ``float``. The type of data to be represented.

    Returns
    -------
    scipy.csr_matrix
        The newly generated matrix
    """
    if data.shape == (0,):
        # an empty vector. Note, this short circuit is necessary as calling
        # csr_matrix([], shape=(0, 0), dtype=dtype) will result in a matrix
        # has a shape of (1, 0).
        return csr_matrix((0, 0), dtype=dtype)
    elif data.shape in ((1, 0), (0, 1)) and data.size == 0:
        # an empty matrix. This short circuit is necessary for the same reason
        # as the empty vector. While a (1, 0) matrix is _empty_, this does
        # confound code that assumes that (1, 0) means there might be metadata
        # or IDs associated with that singular row
        return csr_matrix((0, 0), dtype=dtype)
    elif len(data.shape) == 1:
        # a vector
        shape = (1, data.shape[0])
    else:
        shape = data.shape

    matrix = coo_matrix(data, shape=shape, dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_nparray_to_sparse(data, dtype=float):
    """Takes a list of numpy arrays and creates a scipy.sparse matrix.

    Parameters
    ----------
    data : iterable of numpy.array
        The data to convert into a sparse matrix
    dtype : type, optional
        Defaults to ``float``. The type of data to be represented.

    Returns
    -------
    scipy.csr_matrix
        The newly generated matrix
    """
    matrix = coo_matrix(data, shape=(len(data), len(data[0])), dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_sparse_to_sparse(data, dtype=float):
    """Takes a list of scipy.sparse matrices and creates a scipy.sparse mat.

    Parameters
    ----------
    data : iterable of scipy.sparse matrices
        The data to convert into a sparse matrix
    dtype : type, optional
        Defaults to ``float``. The type of data to be represented.

    Returns
    -------
    scipy.csr_matrix
        The newly generated matrix
    """
    if isspmatrix(data[0]):
        if data[0].shape[0] > data[0].shape[1]:
            n_cols = len(data)
            n_rows = data[0].shape[0]
        else:
            n_rows = len(data)
            n_cols = data[0].shape[1]
    else:
        all_keys = flatten([d.keys() for d in data])
        n_rows = max(all_keys, key=itemgetter(0))[0] + 1
        n_cols = max(all_keys, key=itemgetter(1))[1] + 1
        if n_rows > n_cols:
            n_cols = len(data)
        else:
            n_rows = len(data)

    data = vstack(data)
    matrix = coo_matrix(data, shape=(n_rows, n_cols),
                        dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_dict_to_sparse(data, dtype=float):
    """Takes a list of dict {(row,col):val} and creates a scipy.sparse mat.

    Parameters
    ----------
    data : iterable of dicts
        The data to convert into a sparse matrix
    dtype : type, optional
        Defaults to ``float``. The type of data to be represented.

    Returns
    -------
    scipy.csr_matrix
        The newly generated matrix
    """
    if isspmatrix(data[0]):
        if data[0].shape[0] > data[0].shape[1]:
            is_col = True
            n_cols = len(data)
            n_rows = data[0].shape[0]
        else:
            is_col = False
            n_rows = len(data)
            n_cols = data[0].shape[1]
    else:
        all_keys = flatten([d.keys() for d in data])
        n_rows = max(all_keys, key=itemgetter(0))[0] + 1
        n_cols = max(all_keys, key=itemgetter(1))[1] + 1
        if n_rows > n_cols:
            is_col = True
            n_cols = len(data)
        else:
            is_col = False
            n_rows = len(data)

    rows = []
    cols = []
    vals = []
    for row_idx, row in enumerate(data):
        for (row_val, col_idx), val in row.items():
            if is_col:
                # transpose
                rows.append(row_val)
                cols.append(row_idx)
                vals.append(val)
            else:
                rows.append(row_idx)
                cols.append(col_idx)
                vals.append(val)

    matrix = coo_matrix((vals, (rows, cols)), shape=(n_rows, n_cols),
                        dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def dict_to_sparse(data, dtype=float, shape=None):
    """Takes a dict {(row,col):val} and creates a scipy.sparse matrix.

    Parameters
    ----------
    data : dict
        The data to convert into a sparse matrix
    dtype : type, optional
        Defaults to ``float``. The type of data to be represented.

    Returns
    -------
    scipy.csr_matrix
        The newly generated matrix
    """
    if shape is None:
        n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
        n_cols = max(data.keys(), key=itemgetter(1))[1] + 1
    else:
        n_rows, n_cols = shape

    rows = []
    cols = []
    vals = []
    for (r, c), v in data.items():
        rows.append(r)
        cols.append(c)
        vals.append(v)

    return coo_arrays_to_sparse((vals, (rows, cols)),
                                shape=(n_rows, n_cols), dtype=dtype)
