#!/usr/bin/env python
"""The BIOM Table API"""

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division
import os
import numpy as np
from copy import deepcopy
from datetime import datetime
from json import dumps, loads
from functools import reduce, partial
from operator import itemgetter, add
from itertools import izip
from collections import defaultdict, Hashable
from numpy import ndarray, asarray, zeros, empty, newaxis
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, isspmatrix, vstack

from biom.exception import TableException, UnknownAxisError, UnknownIDError
from biom.util import (get_biom_format_version_string,
                       get_biom_format_url_string, flatten, natsort,
                       prefer_self, index_list, H5PY_VLEN_STR, HAVE_H5PY)

from ._filter import filter_sparse_array
from ._transform import _transform


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Adam Robbins-Pianka",
               "Joshua Shorenstein", "Jose Antonio Navas Molina"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class Table(object):

    """The (canonically pronounced 'teh') Table.

    Give in to the power of the Table!

    """

    def __init__(self, data, observation_ids, sample_ids,
                 observation_metadata=None, sample_metadata=None,
                 table_id=None, type=None, **kwargs):

        self.type = type
        self.table_id = table_id
        self._data = data

        # using object to allow for variable length strings
        self.sample_ids = np.asarray(sample_ids, dtype=object)
        self.observation_ids = np.asarray(observation_ids, dtype=object)

        if sample_metadata is not None:
            self.sample_metadata = tuple(sample_metadata)
        else:
            self.sample_metadata = None

        if observation_metadata is not None:
            self.observation_metadata = tuple(observation_metadata)
        else:
            self.observation_metadata = None

        # These will be set by _index_ids()
        self._sample_index = None
        self._obs_index = None

        self._verify_metadata()
        self._cast_metadata()
        self._index_ids()

    def _index_ids(self):
        """Sets lookups {id:index in _data}.

        Should only be called in constructor as this modifies state.
        """
        self._sample_index = index_list(self.sample_ids)
        self._obs_index = index_list(self.observation_ids)

    def _conv_to_self_type(self, vals, transpose=False, dtype=None):
        """For converting vectors to a compatible self type"""
        if dtype is None:
            dtype = self.dtype

        if isspmatrix(vals):
            return vals
        else:
            return to_sparse(vals, transpose, dtype)

    def _to_dense(self, vec):
        """Converts a row/col vector to a dense numpy array.

        Always returns a 1-D row vector for consistency with numpy iteration
        over arrays.
        """
        dense_vec = np.asarray(vec.todense())

        if vec.shape == (1, 1):
            # Handle the special case where we only have a single element, but
            # we don't want to return a numpy scalar / 0-d array. We still want
            # to return a vector of length 1.
            return dense_vec.reshape(1)
        else:
            return np.squeeze(dense_vec)

    def _verify_metadata(self):
        """Obtain some notion of sanity on object construction with inputs"""
        try:
            n_obs, n_samp = self._data.shape
        except:
            n_obs = n_samp = 0

        if n_obs != len(self.observation_ids):
            raise TableException(
                "Number of observation_ids differs from matrix size!")

        if n_obs != len(set(self.observation_ids)):
            raise TableException("Duplicate observation_ids")

        if n_samp != len(self.sample_ids):
            raise TableException(
                "Number of sample_ids differs from matrix size!")

        if n_samp != len(set(self.sample_ids)):
            raise TableException("Duplicate sample_ids")

        if self.sample_metadata is not None and \
           n_samp != len(self.sample_metadata):
            raise TableException("sample_metadata not in a compatible shape"
                                 "with data matrix!")

        if self.observation_metadata is not None and \
           n_obs != len(self.observation_metadata):
            raise TableException("observation_metadata not in a compatible"
                                 "shape with data matrix!")

    def _cast_metadata(self):
        """Casts all metadata to defaultdict to support default values.

        Should be called after any modifications to sample/observation
        metadata.
        """
        default_samp_md = []
        default_obs_md = []

        # if we have a list of [None], set to None
        if self.sample_metadata is not None:
            if self.sample_metadata.count(None) == len(self.sample_metadata):
                self.sample_metadata = None

        if self.sample_metadata is not None:
            for samp_md in self.sample_metadata:
                d = defaultdict(lambda: None)

                if isinstance(samp_md, dict):
                    d.update(samp_md)
                elif samp_md is None:
                    pass
                else:
                    raise TableException("Unable to cast metadata: %s" %
                                         repr(samp_md))

                default_samp_md.append(d)
            self.sample_metadata = tuple(default_samp_md)

        # if we have a list of [None], set to None
        if self.observation_metadata is not None:
            none_count = self.observation_metadata.count(None)
            if none_count == len(self.observation_metadata):
                self.observation_metadata = None

        if self.observation_metadata is not None:
            for obs_md in self.observation_metadata:
                d = defaultdict(lambda: None)

                if isinstance(obs_md, dict):
                    d.update(obs_md)
                elif obs_md is None:
                    pass
                else:
                    raise TableException("Unable to cast metadata: %s" %
                                         repr(obs_md))

                default_obs_md.append(d)
            self.observation_metadata = tuple(default_obs_md)

    @property
    def shape(self):
        return self._data.shape

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def nnz(self):
        return self._data.nnz

    def add_metadata(self, md, axis='sample'):
        """Take a dict of metadata and add it to an axis.
        Parameters
        ----------
        md : dict of dict
            ``md`` should be of the form ``{id:{dict_of_metadata}}``
        axis : 'sample' or 'observation'
            The axis to operate on
        """
        if axis == 'sample':
            if self.sample_metadata is not None:
                for id_, md_entry in md.iteritems():
                    if self.exists(id_):
                        idx = self.index(id_, 'sample')
                        self.sample_metadata[idx].update(md_entry)
            else:
                self.sample_metadata = tuple([md[id_] if id_ in md else
                                              None for id_ in self.sample_ids])
        elif axis == 'observation':
            if self.observation_metadata is not None:
                for id_, md_entry in md.iteritems():
                    if self.exists(id_, axis="observation"):
                        idx = self.index(id_, 'observation')
                        self.observation_metadata[idx].update(md_entry)
            else:
                self.observation_metadata = tuple([md[id_] if id_ in md else
                                                   None for id_ in
                                                   self.observation_ids])
        else:
            raise UnknownAxisError(axis)

        self._cast_metadata()

    def __getitem__(self, args):
        """Handles row or column slices."""
        if self.is_empty():
            raise IndexError("Cannot retrieve an element from an empty/null "
                             "table.")

        try:
            row, col = args
        except:
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
        """
        self._data = self._data.tocsr()
        return self._data.getrow(row_idx)

    def _get_col(self, col_idx):
        """Return the column at ``col_idx``.

        A column vector will be returned as a scipy.sparse matrix in csc
        format.
        """
        self._data = self._data.tocsc()
        return self._data.getcol(col_idx)

    def reduce(self, f, axis):
        """Reduce over axis with f

        ``axis`` can be either ``sample`` or ``observation``
        """
        if self.is_empty():
            raise TableException("Cannot reduce an empty table")

        # np.apply_along_axis might reduce type conversions here and improve
        # speed. am opting for reduce right now as I think its more readable
        if axis == 'sample':
            return asarray([reduce(f, v) for v in self.iter_data()])
        elif axis == 'observation':
            return asarray([reduce(f, v) for v in
                            self.iter_data(axis="observation")])
        else:
            raise TableException("Unknown reduction axis")

    def sum(self, axis='whole'):
        """Returns the sum by axis

        axis can be:

        ``whole``       : whole matrix sum

        ``sample``     : return a vector with a sum for each sample

        ``observation`` : return a vector with a sum for each observation
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
        """Return a new table that is the transpose of this table.

        The returned table will be an entirely new table, including copies of
        the (transposed) data, sample/observation IDs and metadata.
        """
        sample_md_copy = deepcopy(self.sample_metadata)
        obs_md_copy = deepcopy(self.observation_metadata)

        if self._data.getformat() == 'lil':
            # lil's transpose method doesn't have the copy kwarg, but all of
            # the others do.
            self._data = self._data.tocsr()

        # sample ids and observations are reversed becuase we trasposed
        return self.__class__(self._data.transpose(copy=True),
                              self.sample_ids[:], self.observation_ids[:],
                              sample_md_copy, obs_md_copy, self.table_id)

    def index(self, id_, axis):
        """Return the index of the identified sample/observation.

        Parameters
        ----------
        id_ : str
            ID of the sample or observation whose index will be returned.
        axis : {'sample', 'observation'}
            Axis to search for `id_`.

        Returns
        -------
        int
            Index of the sample/observation identified by `id_`.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        UnknownIDError
            If provided an unrecognized sample/observation ID.

        Examples
        --------

        >>> import numpy as np
        >>> from biom.table import table_factory

        Create a 2x3 BIOM table:

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = table_factory(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Get the index of the observation with ID "O2":

        >>> table.index('O2', 'observation')
        1

        Get the index of the sample with ID "S1":

        >>> table.index('S1', 'sample')
        0

        """
        if axis == 'sample':
            idx_lookup = self._sample_index
        elif axis == 'observation':
            idx_lookup = self._obs_index
        else:
            raise UnknownAxisError(axis)

        if id_ not in idx_lookup:
            raise UnknownIDError(id_, axis)

        return idx_lookup[id_]

    def get_value_by_ids(self, obs_id, samp_id):
        """Return value in the matrix corresponding to ``(obs_id, samp_id)``"""
        return self[self.index(obs_id, 'observation'),
                    self.index(samp_id, 'sample')]

    def __str__(self):
        """Stringify self

        Default str output for a Table is just row/col ids and data values
        """
        return self.delimited_self()

    def exists(self, id_, axis="sample"):
        """Returns whether id_ exists in axis

        Parameters
        ----------
        id_: str
            id to check if exists
        axis : 'sample' or 'observation'
            The axis to check

        Returns
        -------
        bool
            True if ``id_`` exists, False otherwise
        """
        if axis == "sample":
            return id_ in self._sample_index
        elif axis == "observation":
            return id_ in self._obs_index
        else:
            raise UnknownAxisError(axis)

    def delimited_self(self, delim='\t', header_key=None, header_value=None,
                       metadata_formatter=str,
                       observation_column_name='#OTU ID'):
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
        if self.is_empty():
            raise TableException("Cannot delimit self if I don't have data...")

        samp_ids = delim.join(map(str, self.sample_ids))

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
            output = ['# Constructed from biom file',
                      '%s%s%s\t%s' % (observation_column_name, delim, samp_ids,
                                      header_value)]
        else:
            output = ['# Constructed from biom file',
                      '%s%s%s' % (observation_column_name, delim, samp_ids)]

        for obs_id, obs_values in zip(self.observation_ids, self._iter_obs()):
            str_obs_vals = delim.join(map(str, self._to_dense(obs_values)))

            if header_key and self.observation_metadata is not None:
                md = self.observation_metadata[self._obs_index[obs_id]]
                md_out = metadata_formatter(md.get(header_key, None))
                output.append(
                    '%s%s%s\t%s' %
                    (obs_id, delim, str_obs_vals, md_out))
            else:
                output.append('%s%s%s' % (obs_id, delim, str_obs_vals))

        return '\n'.join(output)

    def is_empty(self):
        """Returns ``True`` if the table is empty"""
        if not self.sample_ids.size or not self.observation_ids.size:
            return True
        else:
            return False

    def __iter__(self):
        """Defined by subclass"""
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
        """Returns the fraction of nonzero elements in the table."""
        density = 0.0

        if not self.is_empty():
            density = (self.nnz /
                       (len(self.sample_ids) * len(self.observation_ids)))

        return density

    def descriptive_equality(self, other):
        """For use in testing, describe how the tables are not equal"""
        if self.observation_ids != other.observation_ids:
            return "Observation IDs are not the same"
        if self.sample_ids != other.sample_ids:
            return "Sample IDs are not the same"
        if self.observation_metadata != other.observation_metadata:
            return "Observation metadata are not the same"
        if self.sample_metadata != other.sample_metadata:
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
        if np.any(self.observation_ids != other.observation_ids):
            return False
        if np.any(self.sample_ids != other.sample_ids):
            return False
        if self.observation_metadata != other.observation_metadata:
            return False
        if self.sample_metadata != other.sample_metadata:
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

    def data(self, id_, axis):
        """Returns observations associated with sample id `id_` or
        samples associated with observation id `id_`

        Parameters
        ----------
        id_ : str
            ID of the sample or observation whose index will be returned.
        axis : {'sample', 'observation'}
            Axis to search for `id_`.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        """
        if axis == 'sample':
            return self._to_dense(self[:, self.index(id_, 'sample')])
        elif axis == 'observation':
            return self._to_dense(self[self.index(id_, 'observation'), :])
        else:
            raise UnknownAxisError(axis)

    def copy(self):
        """Returns a copy of the table"""
        # NEEDS TO BE A DEEP COPY, MIGHT NOT GET METADATA! NEED TEST!
        return self.__class__(self._data.copy(),  self.observation_ids[:],
                              self.sample_ids[:], self.observation_metadata,
                              self.sample_metadata, self.table_id)

    def iter_data(self, axis='sample'):
        """Yields axis values

        Parameters
        ----------
        axis: 'sample' or 'observation'
            axis to iterate over

        Returns
        -------
        generator
            yields list of values for each value in 'axis'

        Raises
        ------
        UnknownAxisError
            If axis other than 'sample' or 'observation' passed
        """
        if axis == "sample":
            for samp_v in self._iter_samp():
                yield self._to_dense(samp_v)
        elif axis == "observation":
            for obs_v in self._iter_obs():
                yield self._to_dense(obs_v)
        else:
            raise UnknownAxisError(axis)

    def iter(self, dense=True, axis='sample'):
        """Yields ``(value, id, metadata)``

        NOTE: will return ``None`` in metadata positions if the corresponding
        axis metadata metadata is set to ``None``

        Parameters
        ----------
        dense : bool
            If True, yield values as a dense vector
        axis : str, either 'sample' or 'observation'
            The axis to iterate over

        Returns
        -------
        GeneratorType
            A generator that yields (values, id, metadata)

        """
        if axis == 'sample':
            ids = self.sample_ids
            iter_ = self._iter_samp()
            metadata = self.sample_metadata
        elif axis == 'observation':
            ids = self.observation_ids
            iter_ = self._iter_obs()
            metadata = self.observation_metadata
        else:
            raise UnknownAxisError(axis)

        if metadata is None:
            metadata = (None,) * len(ids)

        for vals, id_, md in izip(iter_, ids, metadata):
            if dense:
                vals = self._to_dense(vals)

            yield (vals, id_, md)

    def sort_order(self, order, axis='sample'):
        """Return a new table with `axis` in `order`

        Parameters
        ----------
        order : iterable
            The desired order for axis
        axis : 'sample' or 'observation'
            The axis to operate on
        """
        md = []
        vals = []
        if axis == 'sample':
            for id_ in order:
                cur_idx = self.index(id_, 'sample')
                vals.append(self._to_dense(self[:, cur_idx]))

                if self.sample_metadata is not None:
                    md.append(self.sample_metadata[cur_idx])

            if not md:
                md = None

            return self.__class__(self._conv_to_self_type(vals,
                                                          transpose=True),
                                  self.observation_ids[:], order[:],
                                  self.observation_metadata, md,
                                  self.table_id)
        elif axis == 'observation':
            for id_ in order:
                cur_idx = self.index(id_, 'observation')
                vals.append(self[cur_idx, :])

                if self.observation_metadata is not None:
                    md.append(self.observation_metadata[cur_idx])

            if not md:
                md = None

            return self.__class__(self._conv_to_self_type(vals),
                                  order[:], self.sample_ids[:],
                                  md, self.sample_metadata, self.table_id)
        else:
            raise UnknownAxisError(axis)

    def sort(self, sort_f=natsort, axis='sample'):
        """Return a table sorted along axis

        Parameters
        ----------
        sort_f : function
            A function that takes a list of values and sorts it
        axis : 'sample' or 'observation'
            The axis to operate on
        """
        if axis == 'sample':
            return self.sort_order(sort_f(self.sample_ids))
        elif axis == 'observation':
            return self.sort_order(sort_f(self.observation_ids),
                                   axis='observation')
        else:
            raise UnknownAxisError(axis)

    def filter(self, ids_to_keep, axis, invert=False):
        """Filter in place a table based on a function or iterable.

        Parameters
        ----------
        ids_to_keep : function or iterable
            If a function, it will be called with the id (a string)
            and the dictionary of metadata of each sample, and must
            return a boolean.
            If it's an iterable, it will be converted to an array of
            bools.
        axis : str
            It controls whether to filter samples or observations. Can
            be "sample" or "observation".
        invert : bool
            If set to True, discard samples or observations where
            `ids_to_keep` returns True
        """
        if axis == 'sample':
            axis = 1
            ids = self.sample_ids
            metadata = self.sample_metadata
        elif axis == 'observation':
            axis = 0
            ids = self.observation_ids
            metadata = self.observation_metadata
        else:
            raise ValueError("Unsupported axis")

        arr = self._data
        arr, ids, metadata = filter_sparse_array(arr,
                                                 ids,
                                                 metadata,
                                                 ids_to_keep,
                                                 axis,
                                                 invert=invert)

        self._data = arr
        if axis == 1:
            self.sample_ids = ids
            self.sample_metadata = metadata
        elif axis == 0:
            self.observation_ids = ids
            self.observation_metadata = metadata

    def partition(self, f, axis='sample'):
        """Yields partitions

        Parameters
        ----------
        f : function
            ``f`` is given the ID and metadata of the vector and must return
            what partition the vector is part of.
        axis : str, either 'sample' or 'observation'
            The axis to iterate over

        Returns
        -------
        GeneratorType
            A generator that yields `Table`

        """
        partitions = {}
        # conversion of vector types is not necessary, vectors are not
        # being passed to an arbitrary function
        for vals, id_, md in self.iter(dense=False, axis=axis):
            part = f(id_, md)

            # try to make it hashable...
            if not isinstance(part, Hashable):
                part = tuple(part)

            if part not in partitions:
                partitions[part] = [[], [], []]

            partitions[part][0].append(id_)
            partitions[part][1].append(vals)
            partitions[part][2].append(md)

        for part, (ids, values, metadata) in partitions.iteritems():
            if axis == 'sample':
                data = self._conv_to_self_type(values, transpose=True)

                samp_ids = ids
                samp_md = metadata
                obs_ids = self.observation_ids[:]

                if self.observation_metadata is not None:
                    obs_md = self.observation_metadata[:]
                else:
                    obs_md = None

            elif axis == 'observation':
                data = self._conv_to_self_type(values, transpose=False)

                obs_ids = ids
                obs_md = metadata
                samp_ids = self.sample_ids[:]

                if self.sample_metadata is not None:
                    samp_md = self.sample_metadata[:]
                else:
                    samp_md = None

            yield part, table_factory(data, obs_ids, samp_ids, obs_md, samp_md,
                                      self.table_id)

    def collapse(self, f, reduce_f=add, norm=True, min_group_size=2,
                 include_collapsed_metadata=True, one_to_many=False,
                 one_to_many_mode='add', one_to_many_md_key='Path',
                 strict=False, axis='sample'):
        """Collapse partitions in a table by metadata or by IDs

        Partition data by metadata or IDs and then collapse each partition into
        a single vector.

        If ``include_collapsed_metadata`` is True, metadata for the collapsed
        partition are retained and can be referred to by the corresponding ID
        from each vector within the partition.

        The remainder is only relevant to setting ``one_to_many`` to True.

        If ``one_to_many`` is True, allow vectors to collapse into multiple
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
        via ``one_to_many_mode``: ``add`` and ``divide``. ``add`` will add the
        vector counts to each partition that the vector maps to, which may
        increase the total number of counts in the output table. ``divide``
        will divide a vectors's counts by the number of metadata that the
        vector has before adding the counts to each partition. This will not
        increase the total number of counts in the output table.

        If ``one_to_many_md_key`` is specified, that becomes the metadata
        key that describes the collapsed path. If a value is not specified,
        then it defaults to 'Path'.

        If ``strict`` is specified, then all metadata pathways operated on
        must be indexable by ``metadata_f``.

        ``one_to_many`` and ``norm`` are not supported together.

        ``one_to_many`` and ``reduce_f`` are not supported together.

        ``one_to_many`` and ``min_group_size`` are not supported together.

        A final note on space consumption. At present, the ``one_to_many``
        functionality requires a temporary dense matrix representation.

        Parameters
        ----------
        f : function
            Function that is used to determine what partition a vector belongs
            to
        reduce_f : function
            Function that reduces two vectors in a one-to-one collapse
        norm : bool
            If `True`, normalize the resulting table
        min_group_size : int
            The minimum size of a partition of performing a one-to-many
            collapse
        include_collapsed_metadata : bool
            If `True`, retain the collapsed metadata keyed by the original IDs
            of the associated vectors
        one_to_many : bool
            Perform a one-to-many collapse
        one_to_many_mode : str, 'add' or 'divide'
            The way to reduce two vectors in a one-to-many collapse
        one_to_many_md_key : str
            The if `include_collapsed_metadata` is `True`, store the original
            vector metadata under this key
        strict : bool
            Requires full pathway data within a one-to-many structure

        Returns
        -------
        Table
            The collapsed table

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
        # axis_slice is only necessary in the one-to-many case
        if axis == 'sample':
            axis_ids_md = lambda t: (t.sample_ids, t.sample_metadata)
            transpose = True
            new_data_shape = lambda ids, collapsed: (len(ids), len(collapsed))
            axis_slice = lambda lookup, key: (slice(None), lookup[key])
        elif axis == 'observation':
            axis_ids_md = lambda t: (t.observation_ids, t.observation_metadata)
            transpose = False
            new_data_shape = lambda ids, collapsed: (len(collapsed), len(ids))
            axis_slice = lambda lookup, key: (lookup[key], slice(None))
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

            for id_, md in izip(*axis_ids_md(self)):
                md_iter = f(id_, md)
                num_md = 0
                while True:
                    try:
                        pathway, partition = md_iter.next()
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
            dtype = float if one_to_many_mode == 'divide' else self.dtype

            new_data = zeros(new_data_shape(axis_ids_md(self)[0], new_md),
                             dtype=dtype)

            # for each vector
            # for each bin in the metadata
            # for each partition associated with the vector
            for vals, id_, md in self.iter(axis=axis):
                md_iter = f(id_, md)

                while True:
                    try:
                        pathway, part = md_iter.next()
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

                    if one_to_many_mode == 'add':
                        new_data[axis_slice(idx_lookup, part)] += vals
                    else:
                        new_data[axis_slice(idx_lookup, part)] += \
                            vals / md_count[id_]

            if include_collapsed_metadata:
                # reassociate pathway information
                for k, i in sorted(idx_lookup.iteritems(), key=itemgetter(1)):
                    collapsed_md.append({one_to_many_md_key: new_md[k]})

            # get the new sample IDs
            collapsed_ids = [k for k, i in sorted(idx_lookup.iteritems(),
                                                  key=itemgetter(1))]

            # convert back to self type
            data = self._conv_to_self_type(new_data)
        else:
            for part, table in self.partition(f, axis=axis):
                axis_ids, axis_md = axis_ids_md(table)

                if len(axis_ids) < min_group_size:
                    continue

                redux_data = table.reduce(reduce_f, self._invert_axis(axis))
                if norm:
                    redux_data /= len(axis_ids)

                collapsed_data.append(self._conv_to_self_type(redux_data))
                collapsed_ids.append(part)

                if include_collapsed_metadata:
                    # retain metadata but store by original id
                    tmp_md = {}
                    for id_, md in izip(axis_ids, axis_md):
                        tmp_md[id_] = md
                    collapsed_md.append(tmp_md)

            data = self._conv_to_self_type(collapsed_data, transpose=transpose)

        # if the table is empty
        if 0 in data.shape:
            raise TableException("Collapsed table is empty!")

        if axis == 'sample':
            sample_ids = collapsed_ids
            sample_md = collapsed_md
            obs_ids = self.observation_ids[:]
            if self.observation_metadata is not None:
                obs_md = self.observation_metadata[:]
            else:
                obs_md = None
        else:
            sample_ids = self.sample_ids[:]
            obs_ids = collapsed_ids
            obs_md = collapsed_md
            if self.sample_metadata is not None:
                sample_md = self.sample_metadata[:]
            else:
                sample_md = None

        return table_factory(data, obs_ids, sample_ids, obs_md, sample_md,
                             self.table_id)

    def _invert_axis(self, axis):
        """Invert an axis"""
        if axis == 'sample':
            return 'observation'
        elif axis == 'observation':
            return 'sample'
        else:
            return UnknownAxisError(axis)

    def transform(self, f, axis='sample'):
        """Iterate over `axis`, applying a function `f` to each vector.

        Only non null values can be modified: the density of the table
        can't increase. However, zeroing values is fine.

        Parameters
        ----------
        f : function
            A function that takes three values: an array of nonzero
            values corresponding to each observation or sample, an
            observation or sample id, and an observation or sample
            metadata entry. It must return an array of transformed
            values that replace the original values.
        axis : 'sample' or 'observation'
            The axis to operate on.
        """
        if axis == 'sample':
            axis = 1
            ids = self.sample_ids
            metadata = self.sample_metadata
            arr = self._data.tocsc()
        elif axis == 'observation':
            axis = 0
            ids = self.observation_ids
            metadata = self.observation_metadata
            arr = self._data.tocsr()
        else:
            raise UnknownAxisError(axis)

        _transform(arr, ids, metadata, f, axis)
        arr.eliminate_zeros()

        self._data = arr

    def norm(self, axis='sample'):
        """Normalize in place sample values by an observation, or vice versa.

        Parameters
        ----------
        axis : 'sample' or 'observation'
            The axis to use for normalization
        """
        def f(val, id_, _):
            return val / float(val.sum())

        self.transform(f, axis=axis)

    def nonzero(self):
        """Returns locations of nonzero elements within the data matrix

        The values returned are ``(observation_id, sample_id)``
        """
        # this is naively implemented. If performance is a concern, private
        # methods can be written to hit against the underlying types directly
        for o_idx, samp_vals in enumerate(self.iter_data(axis="observation")):
            for s_idx in samp_vals.nonzero()[0]:
                yield (self.observation_ids[o_idx], self.sample_ids[s_idx])

    def nonzero_counts(self, axis, binary=False):
        """Get nonzero summaries about an axis

        axis : either 'sample', 'observation', or 'whole'
        binary : sum of nonzero entries, or summing the values of the entries

        Returns a numpy array in index order to the axis
        """
        if binary:
            dtype = 'int'
            op = lambda x: x.nonzero()[0].size
        else:
            dtype = self.dtype
            op = lambda x: x.sum()

        if axis is 'sample':
            # can use np.bincount for CSMat or ScipySparse
            result = zeros(len(self.sample_ids), dtype=dtype)
            for idx, vals in enumerate(self.iter_data()):
                result[idx] = op(vals)
        elif axis is 'observation':
            # can use np.bincount for CSMat or ScipySparse
            result = zeros(len(self.observation_ids), dtype=dtype)
            for idx, vals in enumerate(self.iter_data(axis="observation")):
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

    def merge(self, other, sample='union', observation='union',
              sample_metadata_f=prefer_self,
              observation_metadata_f=prefer_self):
        """Merge two tables together

        The axes, samples and observations, can be controlled independently.
        Both can either work on ``union`` or ``intersection``.

        ``sample_metadata_f`` and ``observation_metadata_f`` define how to
        merge metadata between tables. The default is to just keep the metadata
        associated to self if self has metadata otherwise take metadata from
        other. These functions are given both metadata dictsand must return
        a single metadata dict

        NOTE: There is an implicit type conversion to ``float``. Tables using
        strings as the type are not supported. No check is currently in
        place.

        NOTE: The return type is always that of ``self``
        """
        # determine the sample order in the resulting table
        if sample is 'union':
            new_samp_order = self._union_id_order(self.sample_ids,
                                                  other.sample_ids)
        elif sample is 'intersection':
            new_samp_order = self._intersect_id_order(self.sample_ids,
                                                      other.sample_ids)
        else:
            raise TableException("Unknown sample merge type: %s" % sample)

        # determine the observation order in the resulting table
        if observation is 'union':
            new_obs_order = self._union_id_order(self.observation_ids,
                                                 other.observation_ids)
        elif observation is 'intersection':
            new_obs_order = self._intersect_id_order(self.observation_ids,
                                                     other.observation_ids)
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
        for id_, idx in new_samp_order:
            sample_ids.append(id_)

            # if we have sample metadata, grab it
            if self.sample_metadata is None or not self.exists(id_):
                self_md = None
            else:
                self_md = self.sample_metadata[self_samp_idx[id_]]

            # if we have sample metadata, grab it
            if other.sample_metadata is None or not other.exists(id_):
                other_md = None
            else:
                other_md = other.sample_metadata[other_samp_idx[id_]]

            sample_md.append(sample_metadata_f(self_md, other_md))

        # POSSIBLE DECOMPOSITION
        # resulting observation ids and sample metadata
        obs_ids = []
        obs_md = []
        for id_, idx in new_obs_order:
            obs_ids.append(id_)

            # if we have observation metadata, grab it
            if self.observation_metadata is None or \
               not self.exists(id_, axis="observation"):
                self_md = None
            else:
                self_md = self.observation_metadata[self_obs_idx[id_]]

            # if we have observation metadata, grab it
            if other.observation_metadata is None or \
                    not other.exists(id_, axis="observation"):
                other_md = None
            else:
                other_md = other.observation_metadata[other_obs_idx[id_]]

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
    def from_hdf5(cls, h5grp, order='observation'):
        """Parse an HDF5 formatted BIOM table

        The expected structure of this group is below. A few basic definitions,
        N is the number of observations and M is the number of samples. Data
        are stored in both compressed sparse row (for observation oriented
        operations) and compressed sparse column (for sample oriented
        operations).

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
        order : 'observation' or 'sample' to indicate which data ordering to
            load the table as

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
        if not HAVE_H5PY:
            raise RuntimeError("h5py is not in the environment, HDF5 support "
                               "is not available")

        if order not in ('observation', 'sample'):
            raise ValueError("Unknown order %s!" % order)

        shape = h5grp.attrs['shape']
        type = None if h5grp.attrs['type'] == '' else h5grp.attrs['type']

        # fetch all of the IDs
        obs_ids = h5grp['observation/ids'][:]
        samp_ids = h5grp['sample/ids'][:]

        # fetch all of the metadata
        no_md = np.array(["[]"])
        obs_md = loads(h5grp['observation'].get('metadata', no_md)[0])
        samp_md = loads(h5grp['sample'].get('metadata', no_md)[0])

        # load the data
        data_path = partial(os.path.join, order)
        data = h5grp[data_path("data")]
        indices = h5grp[data_path("indices")]
        indptr = h5grp[data_path("indptr")]
        cs = (data, indices, indptr)

        if order == 'sample':
            matrix = csc_matrix(cs, shape=shape)
        else:
            matrix = csr_matrix(cs, shape=shape)

        return table_factory(matrix, obs_ids, samp_ids,  obs_md or None,
                             samp_md or None, type=type)

    def to_hdf5(self, h5grp, generated_by, compress=True):
        """Store CSC and CSR in place

        The expected structure of this group is below. A few basic definitions,
        N is the number of observations and M is the number of samples. Data
        are stored in both compressed sparse row (for observation oriented
        operations) and compressed sparse column (for sample oriented
        operations).

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
        generated_by : str
        compress : Boolean  'True' means fiels will be compressed with
            gzip, 'False' means no compression

        See Also
        --------
        Table.format_hdf5

        Examples
        --------
        ### is it okay to actually create files in doctest?

        """
        if not HAVE_H5PY:
            raise RuntimeError("h5py is not in the environment, HDF5 support "
                               "is not available")

        def axis_dump(grp, ids, md, order, compression=None):
            """Store for an axis"""
            self._data = self._data.asformat(order)

            len_ids = len(ids)
            len_indptr = len(self._data.indptr)
            len_data = self.nnz

            grp.create_dataset('data', shape=(len_data,),
                               dtype=np.float64,
                               data=self._data.data,
                               compression=compression)
            grp.create_dataset('indices', shape=(len_data,),
                               dtype=np.int32,
                               data=self._data.indices,
                               compression=compression)
            grp.create_dataset('indptr', shape=(len_indptr,),
                               dtype=np.int32,
                               data=self._data.indptr,
                               compression=compression)

            # if we store IDs in the table as numpy arrays then this store
            # is cleaner, as is the parse
            grp.create_dataset('ids', shape=(len_ids,),
                               dtype=H5PY_VLEN_STR,
                               data=[str(i) for i in ids],
                               compression=compression)

            if md is not None:
                md_str = empty(shape=(), dtype=object)
                md_str[()] = dumps(md)
                grp.create_dataset('metadata', shape=(1,),
                                   dtype=H5PY_VLEN_STR,
                                   data=md_str,
                                   compression=compression)

        h5grp.attrs['id'] = self.table_id if self.table_id else "No Table ID"
        h5grp.attrs['type'] = self.type if self.type else ""
        h5grp.attrs['format-url'] = "http://biom-format.org"
        h5grp.attrs['format-version'] = (2, 0)
        h5grp.attrs['generated-by'] = generated_by
        h5grp.attrs['creation-date'] = datetime.now().isoformat()
        h5grp.attrs['shape'] = self.shape
        h5grp.attrs['nnz'] = self.nnz
        compression = None
        if compress is True:
            compression = 'gzip'
        axis_dump(h5grp.create_group('observation'), self.observation_ids,
                  self.observation_metadata, 'csr', compression)
        axis_dump(h5grp.create_group('sample'), self.sample_ids,
                  self.sample_metadata, 'csc', compression)

    def get_biom_format_object(self, generated_by):
        """Returns a dictionary representing the table in BIOM format.

        This dictionary can then be easily converted into a JSON string for
        serialization.

        ``generated_by``: a string describing the software used to build the
        table

        TODO: This method may be very inefficient in terms of memory usage, so
        it needs to be tested with several large tables to determine if
        optimizations are necessary or not (i.e. subclassing JSONEncoder, using
        generators, etc...).
        """
        if (not isinstance(generated_by, str) and
                not isinstance(generated_by, unicode)):
            raise TableException("Must specify a generated_by string")

        # Fill in top-level metadata.
        biom_format_obj = {}
        biom_format_obj["id"] = self.table_id
        biom_format_obj["format"] = get_biom_format_version_string()
        biom_format_obj["format_url"] =\
            get_biom_format_url_string()
        biom_format_obj["generated_by"] = generated_by
        biom_format_obj["date"] = "%s" % datetime.now().isoformat()

        # Determine if we have any data in the matrix, and what the shape of
        # the matrix is.
        try:
            num_rows, num_cols = self.shape
        except:
            num_rows = num_cols = 0
        has_data = True if num_rows > 0 and num_cols > 0 else False

        # Default the matrix element type to test to be an integer in case we
        # don't have any data in the matrix to test.
        test_element = 0
        if has_data:
            test_element = self[0, 0]

        # Determine the type of elements the matrix is storing.
        if isinstance(test_element, int):
            dtype, matrix_element_type = int, "int"
        elif isinstance(test_element, float):
            dtype, matrix_element_type = float, "float"
        elif isinstance(test_element, unicode):
            dtype, matrix_element_type = unicode, "unicode"
        else:
            raise TableException("Unsupported matrix data type.")

        # Fill in details about the matrix.
        biom_format_obj["matrix_element_type"] = "%s" % matrix_element_type
        biom_format_obj["shape"] = [num_rows, num_cols]

        # Fill in details about the rows in the table and fill in the matrix's
        # data.
        biom_format_obj["rows"] = []
        biom_format_obj["data"] = []
        for obs_index, obs in enumerate(self.iter(axis='observation')):
            biom_format_obj["rows"].append(
                {"id": "%s" % obs[1], "metadata": obs[2]})
            # If the matrix is dense, simply convert the numpy array to a list
            # of data values. If the matrix is sparse, we need to store the
            # data in sparse format, as it is given to us in a numpy array in
            # dense format (i.e. includes zeroes) by iter().
            dense_values = list(obs[0])
            sparse_values = []
            for col_index, val in enumerate(dense_values):
                if float(val) != 0.0:
                    sparse_values.append([obs_index, col_index,
                                          dtype(val)])
            biom_format_obj["data"].extend(sparse_values)

        # Fill in details about the columns in the table.
        biom_format_obj["columns"] = []
        for samp in self.iter():
            biom_format_obj["columns"].append(
                {"id": "%s" % samp[1], "metadata": samp[2]})

        return biom_format_obj

    def get_biom_format_json_string(self, generated_by, direct_io=None):
        """Returns a JSON string representing the table in BIOM format.

        ``generated_by``: a string describing the software used to build the
        table

        If direct_io is not None, the final output is written directly to
        direct_io during processing.
        """
        if (not isinstance(generated_by, str) and
                not isinstance(generated_by, unicode)):
            raise TableException("Must specify a generated_by string")

        # Fill in top-level metadata.
        if direct_io:
            direct_io.write('{')
            direct_io.write('"id": "%s",' % str(self.table_id))
            direct_io.write(
                '"format": "%s",' %
                get_biom_format_version_string())
            direct_io.write(
                '"format_url": "%s",' %
                get_biom_format_url_string())
            direct_io.write('"generated_by": "%s",' % generated_by)
            direct_io.write('"date": "%s",' % datetime.now().isoformat())
        else:
            id_ = '"id": "%s",' % str(self.table_id)
            format_ = '"format": "%s",' % get_biom_format_version_string()
            format_url = '"format_url": "%s",' % get_biom_format_url_string()
            generated_by = '"generated_by": "%s",' % generated_by
            date = '"date": "%s",' % datetime.now().isoformat()

        # Determine if we have any data in the matrix, and what the shape of
        # the matrix is.
        try:
            num_rows, num_cols = self.shape
        except:
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
        elif isinstance(test_element, unicode):
            matrix_element_type = "unicode"
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

        # Fill in details about the rows in the table and fill in the matrix's
        # data.
        if direct_io:
            direct_io.write('"data": [')
        else:
            data = ['"data": [']

        max_row_idx = len(self.observation_ids) - 1
        max_col_idx = len(self.sample_ids) - 1
        rows = ['"rows": [']
        have_written = False
        for obs_index, obs in enumerate(self.iter(axis='observation')):
            # i'm crying on the inside
            if obs_index != max_row_idx:
                rows.append('{"id": "%s", "metadata": %s},' % (obs[1],
                                                               dumps(obs[2])))
            else:
                rows.append('{"id": "%s", "metadata": %s}],' % (obs[1],
                                                                dumps(obs[2])))

            # turns out its a pain to figure out when to place commas. the
            # simple work around, at the expense of a little memory
            # (bound by the number of samples) is to build of what will be
            # written, and then add in the commas where necessary.
            built_row = []
            for col_index, val in enumerate(obs[0]):
                if float(val) != 0.0:
                    built_row.append("[%d,%d,%r]" % (obs_index, col_index,
                                                     val))
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
                columns.append('{"id": "%s", "metadata": %s},' % (
                    samp[1], dumps(samp[2])))
            else:
                columns.append('{"id": "%s", "metadata": %s}]' % (
                    samp[1], dumps(samp[2])))

        rows = ''.join(rows)
        columns = ''.join(columns)

        if direct_io:
            direct_io.write(rows)
            direct_io.write(columns)
            direct_io.write('}')
        else:
            return "{%s}" % ''.join([id_, format_, format_url,
                                     generated_by, date,
                                     matrix_element_type, shape,
                                     ''.join(data), rows, columns])

    def get_biom_format_pretty_print(self, generated_by):
        """Returns a 'pretty print' format of a BIOM file

        ``generated_by``: a string describing the software used to build the
        table

        WARNING: This method displays data values in a columnar format and
        can be misleading.
        """
        return dumps(self.get_biom_format_object(generated_by), sort_keys=True,
                     indent=4)


def list_list_to_nparray(data, dtype=float):
    """Convert a list of lists into a nparray

    [[value, value, ..., value], ...]
    """
    return asarray(data, dtype=dtype)


def dict_to_nparray(data, dtype=float):
    """Takes a dict {(row,col):val} and creates a numpy matrix"""
    rows, cols = zip(*data)  # unzip
    mat = zeros((max(rows) + 1, max(cols) + 1), dtype=dtype)

    for (row, col), val in data.iteritems():
        mat[row, col] = val

    return mat


def list_dict_to_nparray(data, dtype=float):
    """Takes a list of dicts {(0,col):val} and creates an numpy matrix

    Expects each dict to represent a row vector
    """
    n_rows = len(data)
    n_cols = max(flatten([d.keys() for d in data]), key=itemgetter(1))[1] + 1

    mat = zeros((n_rows, n_cols), dtype=dtype)

    for row_idx, row in enumerate(data):
        for (_, col_idx), val in row.iteritems():
            mat[row_idx, col_idx] = val

    return mat


def table_factory(data, observation_ids, sample_ids, observation_metadata=None,
                  sample_metadata=None, table_id=None,
                  input_is_dense=False, transpose=False, **kwargs):
    """Construct a table

    Attempts to make 'data' through various means of juggling. Data can be:

        - numpy.array
        - list of numpy.array vectors
        - SparseObj representation
        - dict representation
        - list of SparseObj representation vectors
        - list of lists of sparse values [[row, col, value], ...]
        - list of lists of dense values [[value, value, ...], ...]
        - Scipy COO data (values, (rows, cols))

    Example usage to create a Table object::

        from biom.table import table_factory
        from numpy import array

        sample_ids = ['s1','s2','s3','s4']
        sample_md = [{'pH':4.2,'country':'Peru'},
                     {'pH':5.2,'country':'Peru'},
                     {'pH':5.0,'country':'Peru'},
                     {'pH':4.9,'country':'Peru'}]

        observation_ids = ['o1','o2','o3']
        observation_md = [{'domain':'Archaea'},
                          {'domain':'Bacteria'},
                          {'domain':'Bacteria'}]

        data = array([[1,2,3,4],
                      [-1,6,7,8],
                      [9,10,11,12]])

        t = table_factory(data,
                          observation_ids,
                          sample_ids,
                          observation_md,
                          sample_md)
    """
    if 'dtype' in kwargs:
        dtype = kwargs['dtype']
    else:
        dtype = float

    if 'shape' in kwargs:
        shape = kwargs['shape']
    else:
        shape = None

    # if we have a numpy array
    if isinstance(data, ndarray):
        data = nparray_to_sparse(data, dtype)

    # if we have a list of things
    elif isinstance(data, list):
        if not data:
            raise TableException("No data was supplied. Cannot create "
                                 "an empty table.")

        elif isinstance(data[0], ndarray):
            data = list_nparray_to_sparse(data, dtype)

        elif isinstance(data[0], dict):
            data = list_dict_to_sparse(data, dtype)

        elif isinstance(data[0], list):
            if input_is_dense:
                d = coo_matrix(data)
                data = coo_arrays_to_sparse((d.data, (d.row, d.col)),
                                            dtype=dtype)
            else:
                data = list_list_to_sparse(data, dtype, shape=shape)

        else:
            raise TableException("Unknown nested list type")

    # if we have a dict representation
    elif isinstance(data, dict):
        data = dict_to_sparse(data, dtype)

    elif isinstance(data, tuple) and isinstance(data[0], ndarray):
        data = coo_arrays_to_sparse(data)

    elif isspmatrix(data):
        pass

    else:
        raise TableException("Cannot handle data!")

    return Table(data, observation_ids, sample_ids,
                 sample_metadata=sample_metadata,
                 observation_metadata=observation_metadata,
                 table_id=table_id, **kwargs)


def to_sparse(values, transpose=False, dtype=float):
    """Try to return a populated scipy.sparse matrix.

    NOTE: assumes the max value observed in row and col defines the size of the
    matrix.
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
        mat = dict_to_sparse(values, dtype)
        if transpose:
            mat = mat.T
        return mat
    elif isspmatrix(values):
        mat = values
        if transpose:
            mat = mat.transpose()
        return mat
    else:
        raise TableException("Unknown input type")


def coo_arrays_to_sparse(data, dtype=np.float64, shape=None):
    """Map directly on to the coo_matrix constructor

    data must be (values, (rows, cols))
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

    [[row, col, value], ...]
    """
    rows, cols, values = izip(*data)

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
    """Convert a numpy array to a scipy.sparse matrix."""
    if len(data.shape) == 1:
        shape = (1, data.shape[0])
    else:
        shape = data.shape

    matrix = coo_matrix(data, shape=shape, dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_nparray_to_sparse(data, dtype=float):
    """Takes a list of numpy arrays and creates a scipy.sparse matrix."""
    matrix = coo_matrix(data, shape=(len(data), len(data[0])), dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_sparse_to_sparse(data, dtype=float):
    """Takes a list of scipy.sparse matrices and creates a scipy.sparse mat."""
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

    data = vstack(data)
    matrix = coo_matrix(data, shape=(n_rows, n_cols),
                        dtype=dtype)
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    return matrix


def list_dict_to_sparse(data, dtype=float):
    """Takes a list of dict {(row,col):val} and creates a scipy.sparse mat."""
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


def dict_to_sparse(data, dtype=float):
    """Takes a dict {(row,col):val} and creates a scipy.sparse matrix."""
    n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
    n_cols = max(data.keys(), key=itemgetter(1))[1] + 1

    rows = []
    cols = []
    vals = []
    for (r, c), v in data.iteritems():
        rows.append(r)
        cols.append(c)
        vals.append(v)

    return coo_arrays_to_sparse((vals, (rows, cols)),
                                shape=(n_rows, n_cols), dtype=dtype)
