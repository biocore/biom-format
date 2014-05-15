#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The BIOM Table API"""

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division
import numpy as np
from copy import deepcopy
from datetime import datetime
from json import dumps, loads
from functools import reduce
from operator import itemgetter, add
from itertools import izip
from collections import defaultdict, Hashable
from numpy import ndarray, asarray, zeros, empty, newaxis
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, isspmatrix, vstack

from biom.exception import TableException, UnknownAxisError, UnknownIDError
from biom.util import (get_biom_format_version_string,
                       get_biom_format_url_string, flatten, natsort,
                       prefer_self, index_list, H5PY_VLEN_STR, HAVE_H5PY)

from ._filter import _filter
from ._transform import _transform


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Adam Robbins-Pianka",
               "Joshua Shorenstein", "Jose Antonio Navas Molina",
               "Jorge Cañardo Alastuey"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


MATRIX_ELEMENT_TYPE = {'int': int, 'float': float, 'unicode': unicode,
                       u'int': int, u'float': float, u'unicode': unicode}


class Table(object):

    """The (canonically pronounced 'teh') Table.

    Give in to the power of the Table!

    """

    def __init__(self, data, observation_ids, sample_ids,
                 observation_metadata=None, sample_metadata=None,
                 table_id=None, type=None, create_date=None, generated_by=None,
                 **kwargs):

        self.type = type
        self.table_id = table_id
        self.create_date = create_date
        self.generated_by = generated_by

        if not isspmatrix(data):
            shape = (len(observation_ids), len(sample_ids))
            input_is_dense = kwargs.get('input_is_dense', False)
            self._data = Table._to_sparse(data, input_is_dense=input_is_dense,
                                          shape=shape)
        else:
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
            return Table._to_sparse(vals, transpose, dtype)

    @staticmethod
    def _to_dense(vec):
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
        axis : {'sample', 'observation'}, optional
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

    def metadata(self, id_, axis):
        """Return the metadata of the identified sample/observation.

        Parameters
        ----------
        id_ : str
            ID of the sample or observation whose index will be returned.
        axis : {'sample', 'observation'}
            Axis to search for `id_`.

        Returns
        -------
        defaultdict or None
            The corresponding metadata `defaultdict` or None of that axis does
            not have metadata.

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
            md = self.sample_metadata
        elif axis == 'observation':
            md = self.observation_metadata
        else:
            raise UnknownAxisError(axis)

        idx = self.index(id_, axis=axis)

        return md[idx] if md is not None else None

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

    def exists(self, id_, axis="sample"):
        """Returns whether id_ exists in axis

        Parameters
        ----------
        id_: str
            id to check if exists
        axis : {'sample', 'observation'}, optional
            The axis to check

        Returns
        -------
        bool
            True if ``id_`` exists, False otherwise

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
        return self.__class__(self._data.copy(),
                              self.observation_ids.copy(),
                              self.sample_ids.copy(),
                              deepcopy(self.observation_metadata),
                              deepcopy(self.sample_metadata),
                              self.table_id)

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
        axis : {'sample', 'observation'}, optional
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
        axis : {'sample', 'observation'}, optional
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
        axis : {'sample', 'observation'}, optional
            The axis to operate on
        """
        if axis == 'sample':
            return self.sort_order(sort_f(self.sample_ids))
        elif axis == 'observation':
            return self.sort_order(sort_f(self.observation_ids),
                                   axis='observation')
        else:
            raise UnknownAxisError(axis)

    def filter(self, ids_to_keep, axis='sample', invert=False, inplace=True):
        """Filter a table based on a function or iterable.

        Parameters
        ----------
        ids_to_keep : function(id, metadata, values) -> bool, or iterable
            If a function, it will be called with the id (a string),
            the dictionary of metadata of each sample/observation and
            the nonzero values of the sample/observation, and must
            return a boolean.
            If it's an iterable, it will be converted to an array of
            bools.
        axis : {'sample', 'observation'}, optional
            It controls whether to filter samples or observations. Can
            be "sample" or "observation".
        invert : bool
            If set to True, discard samples or observations where
            `ids_to_keep` returns True
        inplace : bool, defaults to True
            Whether to return a new table or modify itself.

        Returns
        -------
        biom.Table
            Returns itself if `inplace`, else returns a new filtered table.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.

        """
        table = self if inplace else self.copy()

        if axis == 'sample':
            axis = 1
            ids = table.sample_ids
            metadata = table.sample_metadata
        elif axis == 'observation':
            axis = 0
            ids = table.observation_ids
            metadata = table.observation_metadata
        else:
            raise UnknownAxisError(axis)

        arr = table._data
        arr, ids, metadata = _filter(arr,
                                     ids,
                                     metadata,
                                     ids_to_keep,
                                     axis,
                                     invert=invert)

        table._data = arr
        if axis == 1:
            table.sample_ids = ids
            table.sample_metadata = metadata
        elif axis == 0:
            table.observation_ids = ids
            table.observation_metadata = metadata

        return table

    def partition(self, f, axis='sample'):
        """Yields partitions

        Parameters
        ----------
        f : function
            ``f`` is given the ID and metadata of the vector and must return
            what partition the vector is part of.
        axis : {'sample', 'observation'}, optional
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

            yield part, Table(data, obs_ids, samp_ids, obs_md, samp_md,
                              self.table_id, type=self.type)

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

    def pa(self, inplace=True):
        """Convert the `Table` to presence/absence data

        Parameters
        ----------
        inplace : bool, optional
            Defaults to `False`

        Returns
        -------
        biom.Table
            Returns itself if `inplace`, else returns a new presence/absence
            table.

        Examples
        --------
        >>> import numpy as np
        >>> from biom.table import Table

        Create a 2x3 BIOM table

        >>> data = np.asarray([[0, 0, 1], [1, 3, 42]])
        >>> table = table = Table(data, ['O1', 'O2'], ['S1', 'S2', 'S3'])

        Convert to presence/absence data

        >>> _ = table.pa()
        >>> print table.data('O1', 'observation')
        [ 0.  0.  1.]
        >>> print table.data('O2', 'observation')
        [ 1.  1.  1.]

        """
        def transform_f(data, id_, metadata):
            return np.ones(len(data), dtype=float)

        return self.transform(transform_f, inplace=inplace)

    def transform(self, f, axis='sample', inplace=True):
        """Iterate over `axis`, applying a function `f` to each vector.

        Only non null values can be modified: the density of the table
        can't increase. However, zeroing values is fine.

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
        inplace : bool, defaults to True
            Whether to return a new table or modify itself.

        Returns
        -------
        biom.Table
            Returns itself if `inplace`, else returns a new transformed table.

        Raises
        ------
        UnknownAxisError
            If provided an unrecognized axis.
        """
        table = self if inplace else self.copy()

        if axis == 'sample':
            axis = 1
            ids = table.sample_ids
            metadata = table.sample_metadata
            arr = table._data.tocsc()
        elif axis == 'observation':
            axis = 0
            ids = table.observation_ids
            metadata = table.observation_metadata
            arr = table._data.tocsr()
        else:
            raise UnknownAxisError(axis)

        _transform(arr, ids, metadata, f, axis)
        arr.eliminate_zeros()

        table._data = arr

        return table

    def norm(self, axis='sample'):
        """Normalize in place sample values by an observation, or vice versa.

        Parameters
        ----------
        axis : {'sample', 'observation'}, optional
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
        ./id                         : str, an arbitrary ID
        ./type                       : str, the table type (e.g, OTU table)
        ./format-url                 : str, a URL that describes the format
        ./format-version             : two element tuple of int32,
                                       major and minor
        ./generated-by               : str, what generated this file
        ./creation-date              : str, ISO format
        ./shape                      : two element tuple of int32, N by M
        ./nnz                        : int32 or int64, number of non zero elems
        ./observation                : Group
        ./observation/ids            : (N,) dataset of str or vlen str
        ./observation/matrix         : Group
        ./observation/matrix/data    : (N,) dataset of float64
        ./observation/matrix/indices : (N,) dataset of int32
        ./observation/matrix/indptr  : (M+1,) dataset of int32
        [./observation/metadata]     : Optional, JSON str, in index order
                                       with ids
        ./sample                     : Group
        ./sample/ids                 : (M,) dataset of str or vlen str
        ./sample/matrix              : Group
        ./sample/matrix/data         : (M,) dataset of float64
        ./sample/matrix/indices      : (M,) dataset of int32
        ./sample/matrix/indptr       : (N+1,) dataset of int32
        [./sample/metadata]          : Optional, JSON str, in index order
                                       with ids

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

        create_date = h5grp.attrs['creation-date']
        generated_by = h5grp.attrs['generated-by']

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
        data_grp = h5grp[order]['matrix']
        data = data_grp["data"]
        indices = data_grp["indices"]
        indptr = data_grp["indptr"]
        cs = (data, indices, indptr)

        if order == 'sample':
            matrix = csc_matrix(cs, shape=shape)
        else:
            matrix = csr_matrix(cs, shape=shape)

        return Table(matrix, obs_ids, samp_ids,  obs_md or None,
                     samp_md or None, type=type, create_date=create_date,
                     generated_by=generated_by)

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
        ./id                         : str, an arbitrary ID
        ./type                       : str, the table type (e.g, OTU table)
        ./format-url                 : str, a URL that describes the format
        ./format-version             : two element tuple of int32,
                                       major and minor
        ./generated-by               : str, what generated this file
        ./creation-date              : str, ISO format
        ./shape                      : two element tuple of int32, N by M
        ./nnz                        : int32 or int64, number of non zero elems
        ./observation                : Group
        ./observation/ids            : (N,) dataset of str or vlen str
        ./observation/matrix         : Group
        ./observation/matrix/data    : (N,) dataset of float64
        ./observation/matrix/indices : (N,) dataset of int32
        ./observation/matrix/indptr  : (M+1,) dataset of int32
        [./observation/metadata]     : Optional, JSON str, in index order
                                       with ids
        ./sample                     : Group
        ./sample/ids                 : (M,) dataset of str or vlen str
        ./sample/matrix              : Group
        ./sample/matrix/data         : (M,) dataset of float64
        ./sample/matrix/indices      : (M,) dataset of int32
        ./sample/matrix/indptr       : (N+1,) dataset of int32
        [./sample/metadata]          : Optional, JSON str, in index order
                                       with ids

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

    @classmethod
    def from_json(self, json_table, data_pump=None,
                  input_is_dense=False):
        """Parse a biom otu table type"""
        sample_ids = [col['id'] for col in json_table['columns']]
        sample_metadata = [col['metadata'] for col in json_table['columns']]
        obs_ids = [row['id'] for row in json_table['rows']]
        obs_metadata = [row['metadata'] for row in json_table['rows']]
        dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

        if data_pump is None:
            table_obj = Table(json_table['data'], obs_ids, sample_ids,
                              obs_metadata, sample_metadata,
                              shape=json_table['shape'],
                              dtype=dtype,
                              input_is_dense=input_is_dense)
        else:
            table_obj = Table(data_pump, obs_ids, sample_ids,
                              obs_metadata, sample_metadata,
                              shape=json_table['shape'],
                              dtype=dtype,
                              input_is_dense=input_is_dense)

        return table_obj

    def to_json(self, generated_by, direct_io=None):
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
                rows.append('{"id": %s, "metadata": %s},' % (dumps(obs[1]),
                                                             dumps(obs[2])))
            else:
                rows.append('{"id": %s, "metadata": %s}],' % (dumps(obs[1]),
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
                columns.append('{"id": %s, "metadata": %s},' % (
                    dumps(samp[1]), dumps(samp[2])))
            else:
                columns.append('{"id": %s, "metadata": %s}]' % (
                    dumps(samp[1]), dumps(samp[2])))

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

    @staticmethod
    def from_tsv(lines, obs_mapping, sample_mapping,
                 process_func, **kwargs):
        """Parse an tab separated (observation x sample) formatted BIOM table

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
        Table
            A BIOM ``Table`` object

        Examples
        --------
        >>> from biom.table import Table
        >>> from StringIO import StringIO
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
    def _extract_data_from_tsv(lines, delim='\t', dtype=float,
                               header_mark=None, md_parse=None):
        """Parse a classic table into (sample_ids, obs_ids, data, metadata,
                                       name)
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

        Parameters
        ----------
        lines: list or file-like object
            delimted data to parse
        delim: string
            delimeter in file lines
        dtype: type
        header_mark:  string or None
            string that indicates start of header line
        md_parse:  function or None
            funtion used to parse metdata

        Notes
        ------
        This is intended to be close to how QIIME classic OTU tables are
        parsed with the exception of the additional md_name field

        This function is ported from QIIME (http://www.qiime.org), previously
        named
        parse_classic_otu_table. QIIME is a GPL project, but we obtained
        permission
        from the authors of this function to port it to the BIOM Format project
        (and keep it under BIOM's BSD license).
        """
        if not isinstance(lines, list):
            try:
                lines = lines.readlines()
            except AttributeError:
                raise RuntimeError(
                    "Input needs to support readlines or be indexable")

        # find header, the first line that is not empty and does not start
        # with a #
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

    def to_tsv(self, header_key=None, header_value=None,
               metadata_formatter=str, observation_column_name='#OTU ID'):
        """Return self as a string in tab delimited form

        Default ``str`` output for the ``Table`` is just row/col ids and table
        data without any metadata

        Parameters
        ----------
        header_key : string or None
        header_value : string or None
        metadata_formatter : function
            a function which takes a metadata entry and
            returns a formatted version that should be written to file
        observation_column_name : string
            the name of the first column in the output
            table, corresponding to the observation IDs.

        Returns
        -------
        string
            tab delimited represtation of the Table
            For example, the default will look something like:
                #OTU ID\tSample1\tSample2
                OTU1\t10\t2
                OTU2\t4\t8

        """
        return self.delimited_self('\t', header_key, header_value,
                                   metadata_formatter,
                                   observation_column_name)


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


def dict_to_sparse(data, dtype=float, shape=None):
    """Takes a dict {(row,col):val} and creates a scipy.sparse matrix."""
    if shape is None:
        n_rows = max(data.keys(), key=itemgetter(0))[0] + 1
        n_cols = max(data.keys(), key=itemgetter(1))[1] + 1
    else:
        n_rows, n_cols = shape

    rows = []
    cols = []
    vals = []
    for (r, c), v in data.iteritems():
        rows.append(r)
        cols.append(c)
        vals.append(v)

    return coo_arrays_to_sparse((vals, (rows, cols)),
                                shape=(n_rows, n_cols), dtype=dtype)
