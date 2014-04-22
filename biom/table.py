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

import numpy as np
from copy import deepcopy
from datetime import datetime
from json import dumps
from operator import itemgetter, xor, add
from itertools import izip
from collections import defaultdict, Hashable
from numpy import ndarray, asarray, zeros, empty
import h5py

from biom import get_sparse_backend
from biom.exception import TableException, UnknownID
from biom.util import (get_biom_format_version_string,
                       get_biom_format_url_string, flatten, natsort,
                       prefer_self, index_list)
from functools import reduce

# Define a variable length string type
H5PY_VLEN_STR = h5py.special_dtype(vlen=str)

SparseObj, to_sparse, dict_to_sparseobj, list_dict_to_sparseobj, \
    list_nparray_to_sparseobj, nparray_to_sparseobj, \
    list_list_to_sparseobj = get_sparse_backend()

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Adam Robbins-Pianka"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class Table(object):

    """Abstract base class representing a Table.

    Once created, a Table object is immutable except for its sample/observation
    metadata, which can be modified in place via add_sample_metadata and
    add_observation_metadata.

    Code to simulate immutability taken from here:
        http://en.wikipedia.org/wiki/Immutable_object
    """

    def __setattr__(self, *args):
        raise TypeError("A Table object cannot be modified once created.")
    __delattr__ = __setattr__

    def __init__(self, data, sample_ids, observation_ids, sample_metadata=None,
                 observation_metadata=None, table_id=None,
                 type=None, **kwargs):
        if type is None:
            type = 'Unspecified'

        super(Table, self).__setattr__('type', type)
        super(Table, self).__setattr__('table_id', table_id)
        super(Table, self).__setattr__('_data', data)
        super(Table, self).__setattr__('_dtype', data.dtype)

        # Cast to tuple for immutability.
        super(Table, self).__setattr__('sample_ids', tuple(sample_ids))
        super(Table, self).__setattr__('observation_ids',
                                       tuple(observation_ids))

        if sample_metadata is not None:
            super(Table, self).__setattr__('sample_metadata',
                                           tuple(sample_metadata))
        else:
            super(Table, self).__setattr__('sample_metadata', None)

        if observation_metadata is not None:
            super(Table, self).__setattr__('observation_metadata',
                                           tuple(observation_metadata))
        else:
            super(Table, self).__setattr__('observation_metadata', None)

        # These will be set by _index_ids()
        super(Table, self).__setattr__('_sample_index', None)
        super(Table, self).__setattr__('_obs_index', None)

        self._verify_metadata()
        self._cast_metadata()
        self._index_ids()

    def _index_ids(self):
        """Sets lookups {id:index in _data}.

        Should only be called in constructor as this modifies state.
        """
        super(Table, self).__setattr__('_sample_index',
                                       index_list(self.sample_ids))
        super(Table, self).__setattr__('_obs_index',
                                       index_list(self.observation_ids))

    def _conv_to_self_type(self, vals, transpose=False, dtype=None):
        """For converting vectors to a compatible self type"""
        if dtype is None:
            dtype = self._dtype

        if isinstance(vals, self._data.__class__):
            return vals
        else:
            return to_sparse(vals, transpose, dtype)

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
            raise TableException("sample_metadata not in a compatible shape \
                                   with data matrix!")

        if self.observation_metadata is not None and \
           n_obs != len(self.observation_metadata):
            raise TableException("observation_metadata not in a compatible \
                                   shape with data matrix!")

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
                super(Table, self).__setattr__('sample_metadata', None)

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
            super(Table, self).__setattr__('sample_metadata',
                                           tuple(default_samp_md))

        # if we have a list of [None], set to None
        if self.observation_metadata is not None:
            none_count = self.observation_metadata.count(None)
            if none_count == len(self.observation_metadata):
                super(Table, self).__setattr__('observation_metadata', None)

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
            super(Table, self).__setattr__('observation_metadata',
                                           tuple(default_obs_md))

    def add_observation_metadata(self, md):
        """Take a dict of metadata and add it to an observation.

        ``md`` should be of the form ``{observation_id:{dict_of_metadata}}``
        """
        if self.observation_metadata is not None:
            for id_, md_entry in md.items():
                if self.observation_exists(id_):
                    self.observation_metadata[
                        self.get_observation_index(id_)].update(md_entry)
        else:
            super(Table, self).__setattr__('observation_metadata',
                                           tuple([md[id_] if id_ in md else
                                                  None for id_ in
                                                  self.observation_ids]))
        self._cast_metadata()

    def add_sample_metadata(self, md):
        """Take a dict of metadata and add it to a sample.

        ``md`` should be of the form ``{sample_id:{dict_of_metadata}}``
        """
        if self.sample_metadata is not None:
            for id_, md_entry in md.items():
                if self.sample_exists(id_):
                    self.sample_metadata[
                        self.get_sample_index(id_)].update(md_entry)
        else:
            super(Table, self).__setattr__('sample_metadata',
                                           tuple([md[id_] if id_ in md else
                                                  None
                                                  for id_ in self.sample_ids]))
        self._cast_metadata()

    def __getitem__(self, args):
        """Passes through to internal matrix"""
        return self._data[args]

    def reduce(self, f, axis):
        """Reduce over axis with f

        ``axis`` can be either ``sample`` or ``observation``
        """
        if self.is_empty():
            raise TableException("Cannot reduce an empty table")

        # np.apply_along_axis might reduce type conversions here and improve
        # speed. am opting for reduce right now as I think its more readable
        if axis == 'sample':
            return asarray([reduce(f, v) for v in self.iter_sample_data()])
        elif axis == 'observation':
            return asarray([reduce(f, v) for v in
                            self.iter_observation_data()])
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
            return self._data.sum()
        elif axis == 'sample':
            return self._data.sum(axis=0)
        elif axis == 'observation':
            return self._data.sum(axis=1)
        else:
            raise TableException("Unknown axis '%s'" % axis)

    def transpose(self):
        """Return a new table that is the transpose of this table.

        The returned table will be an entirely new table, including copies of
        the (transposed) data, sample/observation IDs and metadata.
        """
        sample_md_copy = deepcopy(self.sample_metadata)
        obs_md_copy = deepcopy(self.observation_metadata)

        return self.__class__(self._data.T, self.observation_ids[:],
                              self.sample_ids[:], obs_md_copy, sample_md_copy,
                              self.table_id)

    def get_sample_index(self, samp_id):
        """Returns the sample index for sample ``samp_id``"""
        if samp_id not in self._sample_index:
            raise UnknownID("SampleId %s not found!" % samp_id)
        return self._sample_index[samp_id]

    def get_observation_index(self, obs_id):
        """Returns the observation index for observation ``obs_id``"""
        if obs_id not in self._obs_index:
            raise UnknownID("ObservationId %s not found!" % obs_id)
        return self._obs_index[obs_id]

    def get_value_by_ids(self, obs_id, samp_id):
        """Return value in the matrix corresponding to ``(obs_id, samp_id)``
        """
        if obs_id not in self._obs_index:
            raise UnknownID("ObservationId %s not found!" % obs_id)
        if samp_id not in self._sample_index:
            raise UnknownID("SampleId %s not found!" % samp_id)

        return self._data[self._obs_index[obs_id], self._sample_index[samp_id]]

    def __str__(self):
        """Stringify self

        Default str output for a Table is just row/col ids and data values
        """
        return self.delimited_self()

    def sample_exists(self, id_):
        """Returns True if sample ``id_`` exists, False otherwise"""
        return id_ in self._sample_index

    def observation_exists(self, id_):
        """Returns True if observation ``id_`` exists, False otherwise"""
        return id_ in self._obs_index

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
            str_obs_vals = delim.join(map(str, self._conv_to_np(obs_values)))

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
        if not self.sample_ids or not self.observation_ids:
            return True
        else:
            return False

    def __iter__(self):
        """Defined by subclass"""
        return self.iter_samples()

    def _iter_samp(self):
        """Return sample vectors of data matrix vectors"""
        rows, cols = self._data.shape
        for c in range(cols):
            # this pulls out col vectors but need to convert to the expected
            # row vector
            colvec = self._data.get_col(c)
            yield colvec.T

    def _iter_obs(self):
        """Return observation vectors of data matrix"""
        for r in range(self._data.shape[0]):
            yield self._data.get_row(r)

    def get_table_density(self):
        """Returns the fraction of nonzero elements in the table."""
        density = 0.0

        if not self.is_empty():
            density = (self._data.size / (len(self.sample_ids) *
                                          len(self.observation_ids)))

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
        if not self._data_equality(other):
            return "Data elements are not the same"

        return "Tables appear equal"

    def __eq__(self, other):
        """Equality is determined by the data matrix, metadata, and IDs"""
        if self.observation_ids != other.observation_ids:
            return False
        if self.sample_ids != other.sample_ids:
            return False
        if self.observation_metadata != other.observation_metadata:
            return False
        if self.sample_metadata != other.sample_metadata:
            return False
        if not self._data_equality(other):
            return False

        return True

    def _data_equality(self, other):
        """Two SparseObj matrices are equal if the items are equal"""
        if isinstance(self, other.__class__):
            return sorted(self._data.items()) == sorted(other._data.items())

        for s_v, o_v in izip(self.iter_sample_data(),
                             other.iter_sample_data()):
            if not (s_v == o_v).all():
                return False

        return True

    def __ne__(self, other):
        return not (self == other)

    def _conv_to_np(self, v):
        """Converts a vector to a numpy array

        Always returns a row vector for consistancy with numpy iteration over
        arrays
        """
        return SparseObj.convert_vector_to_dense(v)

    def sample_data(self, id_):
        """Return observations associated with sample id ``id_``"""
        if id_ not in self._sample_index:
            raise UnknownID("ID %s is not a known sample ID!" % id_)
        return self._conv_to_np(self._data[:, self._sample_index[id_]])

    def observation_data(self, id_):
        """Return samples associated with observation id ``id_``"""
        if id_ not in self._obs_index:
            raise UnknownID("ID %s is not a known observation ID!" % id_)
        return self._conv_to_np(self._data[self._obs_index[id_], :])

    def copy(self):
        """Returns a copy of the table"""
        # NEEDS TO BE A DEEP COPY, MIGHT NOT GET METADATA! NEED TEST!
        return self.__class__(self._data.copy(), self.sample_ids[:],
                              self.observation_ids[:], self.sample_metadata,
                              self.observation_metadata, self.table_id)

    def iter_sample_data(self):
        """Yields sample values"""
        for samp_v in self._iter_samp():
            yield self._conv_to_np(samp_v)

    def iter_observation_data(self):
        """Yields observation values"""
        for obs_v in self._iter_obs():
            yield self._conv_to_np(obs_v)

    def iter_samples(self, conv_to_np=True):
        """Yields ``(sample_value, sample_id, sample_metadata)``

        NOTE: will return ``None`` in ``sample_metadata`` positions if
        ``self.sample_metadata`` is set to ``None``
        """
        if self.sample_metadata is None:
            samp_metadata = (None,) * len(self.sample_ids)
        else:
            samp_metadata = self.sample_metadata

        iterator = izip(self._iter_samp(), self.sample_ids, samp_metadata)
        for samp_v, samp_id, samp_md in iterator:
            if conv_to_np:
                yield (self._conv_to_np(samp_v), samp_id, samp_md)
            else:
                yield (samp_v, samp_id, samp_md)

    def iter_observations(self, conv_to_np=True):
        """Yields ``(observation_value, observation_id, observation_metadata)``

        NOTE: will return ``None`` in ``observation_metadata`` positions if
        ``self.observation_metadata`` is set to ``None``
        """
        if self.observation_metadata is None:
            obs_metadata = (None,) * len(self.observation_ids)
        else:
            obs_metadata = self.observation_metadata

        iterator = izip(self._iter_obs(), self.observation_ids, obs_metadata)
        for obs_v, obs_id, obs_md in iterator:
            if conv_to_np:
                yield (self._conv_to_np(obs_v), obs_id, obs_md)
            else:
                yield (obs_v, obs_id, obs_md)

    def sort_sample_order(self, sample_order):
        """Return a new table with samples in ``sample_order``"""
        samp_md = []
        vals = []

        for id_ in sample_order:
            cur_idx = self._sample_index[id_]
            vals.append(self._conv_to_np(self[:, cur_idx]))

            if self.sample_metadata is not None:
                samp_md.append(self.sample_metadata[cur_idx])

        if not samp_md:
            samp_md = None

        return self.__class__(self._conv_to_self_type(vals, transpose=True),
                              sample_order[:], self.observation_ids[:],
                              samp_md,
                              self.observation_metadata, self.table_id)

    def sort_observation_order(self, obs_order):
        """Return a new table with observations in ``observation order``"""
        obs_md = []
        vals = []

        for id_ in obs_order:
            cur_idx = self._obs_index[id_]
            vals.append(self[cur_idx, :])

            if self.observation_metadata is not None:
                obs_md.append(self.observation_metadata[cur_idx])

        if not obs_md:
            obs_md = None

        return self.__class__(self._conv_to_self_type(vals),
                              self.sample_ids[:], obs_order[
                                  :], self.sample_metadata,
                              obs_md, self.table_id)

    def sort_by_sample_id(self, sort_f=natsort):
        """Return a table where samples are sorted by ``sort_f``

            ``sort_f`` must take a single parameter: the list of sample ids
        """
        return self.sort_sample_order(sort_f(self.sample_ids))

    def sort_by_observation_id(self, sort_f=natsort):
        """Return a table where observations are sorted by ``sort_f``

            ``sort_f`` must take a single parameter: the list of observation
            ids
        """
        return self.sort_observation_order(sort_f(self.observation_ids))

    # a good refactor in the future is a general filter() method and then
    # specify the axis, like Table.reduce

    # take() is tempting here as well...
    def filter_samples(self, f, invert=False):
        """Filter samples from self based on ``f``

        ``f`` must accept three variables, the sample values, sample ID and
        sample metadata. The function must only return ``True`` or ``False``,
        where ``True`` indicates that a sample should be retained.

        invert: if ``invert == True``, a return value of ``True`` from ``f``
        indicates that a sample should be discarded
        """
        samp_ids = []
        samp_vals = []
        samp_metadata = []

        # builtin filter puts all of this into memory and then return to the
        # for loop. This will impact memory substantially on large sparse
        # matrices
        for s_val, s_id, s_md in self.iter_samples():
            if not xor(f(s_val, s_id, s_md), invert):
                continue

            # there is an implicit converstion to numpy types, want to make
            # sure to convert back to underlying representation.
            samp_vals.append(self._conv_to_self_type(s_val))
            samp_metadata.append(s_md)
            samp_ids.append(s_id)

        # if we don't have any values to keep, throw an exception as we can
        # create an inconsistancy in which there are observation ids but no
        # matrix data in the resulting table
        if not samp_ids:
            raise TableException("All samples were filtered out!")

        # the additional call to _conv_to_self_type is to convert a list of
        # vectors to a matrix
        # transpose is necessary as the underlying storage is sample == col
        return self.__class__(
            self._conv_to_self_type(samp_vals, transpose=True),
            samp_ids[:], self.observation_ids[
                :], samp_metadata,
            self.observation_metadata, self.table_id)

    def filter_observations(self, f, invert=False):
        """Filter observations from self based on ``f``

        ``f`` must accept three variables, the observation values, observation
        ID and observation metadata. The function must only return ``True`` or
        ``False``, where ``True`` indicates that an observation should be
        retained.

        invert: if ``invert == True``, a return value of ``True`` from ``f``
        indicates that an observation should be discarded
        """
        obs_ids = []
        obs_vals = []
        obs_metadata = []

        # builtin filter puts all of this into memory and then return to the
        # for loop. This will impact memory substantially on large sparse
        # matrices
        for o_val, o_id, o_md in self.iter_observations():
            if not xor(f(o_val, o_id, o_md), invert):
                continue

            # there is an implicit converstion to numpy types, want to make
            # sure to convert back to underlying representation.
            obs_vals.append(self._conv_to_self_type(o_val))
            obs_metadata.append(o_md)
            obs_ids.append(o_id)

        # if we don't have any values to keep, throw an exception as we can
        # create an inconsistancy in which there are sample ids but no
        # matrix data in the resulting table
        if not obs_vals:
            raise TableException("All observations were filtered out!")

        return self.__class__(
            self._conv_to_self_type(obs_vals), self.sample_ids[:],
            obs_ids[:], self.sample_metadata, obs_metadata, self.table_id)

    def bin_samples_by_metadata(self, f, constructor=None):
        """Yields tables by metadata

        ``f`` is given the sample metadata by row and must return what "bin"
        the sample is part of.

        ``constructor``: the type of binned tables to create, e.g.
        SparseTaxonTable. If None, the binned tables will be the same type as
        this table.
        """
        if constructor is None:
            constructor = self.__class__

        bins = {}
        # conversion of vector types is not necessary, vectors are not
        # being passed to an arbitrary function
        for samp_v, samp_id, samp_md in self.iter_samples(conv_to_np=False):
            bin = f(samp_md)

            # try to make it hashable...
            if not isinstance(bin, Hashable):
                bin = tuple(bin)

            if bin not in bins:
                bins[bin] = [[], [], []]

            bins[bin][0].append(samp_id)
            bins[bin][1].append(samp_v)
            bins[bin][2].append(samp_md)

        for bin, (samp_ids, samp_values, samp_md) in bins.iteritems():
            data = self._conv_to_self_type(samp_values, transpose=True)
            yield bin, table_factory(data, samp_ids[:],
                                     self.observation_ids[:],
                                     samp_md, self.observation_metadata,
                                     self.table_id,
                                     constructor=constructor)

    def bin_observations_by_metadata(self, f, constructor=None):
        """Yields tables by metadata

        ``f`` is given the observation metadata by row and must return what
        "bin" the observation is part of.

        ``constructor``: the type of binned tables to create, e.g.
        SparseTaxonTable. If None, the binned tables will be the same type as
        this table.
        """
        if constructor is None:
            constructor = self.__class__

        bins = {}
        # conversion of vector types is not necessary, vectors are not
        # being passed to an arbitrary function
        for obs_v, obs_id, obs_md in self.iter_observations(conv_to_np=False):
            bin = f(obs_md)

            # try to make it hashable...
            if not isinstance(bin, Hashable):
                bin = tuple(bin)

            if bin not in bins:
                bins[bin] = [[], [], []]

            bins[bin][0].append(obs_id)
            bins[bin][1].append(obs_v)
            bins[bin][2].append(obs_md)

        for bin, (obs_ids, obs_values, obs_md) in bins.iteritems():
            yield bin, table_factory(self._conv_to_self_type(obs_values),
                                     self.sample_ids[:], obs_ids[
                                         :], self.sample_metadata,
                                     obs_md, self.table_id,
                                     constructor=constructor)

    def collase_samples_by_metadata(self, metadata_f, reduce_f=add, norm=True,
                                    min_group_size=2,
                                    include_collapsed_metadata=True,
                                    constructor=None, one_to_many=False,
                                    one_to_many_mode='add',
                                    one_to_many_md_key='Path', strict=False):
        """Collapse samples in a table by sample metadata

        Bin samples by metadata then collapse each bin into a single sample.

        If ``include_collapsed_metadata`` is True, metadata for the collapsed
        samples are retained and can be referred to by the ``SampleId`` from
        each sample within the bin.

        ``constructor``: the type of collapsed table to create, e.g.
        SparseTaxonTable. If None, the collapsed table will be the same type as
        this table.

        The remainder is only relevant to setting ``one_to_many`` to True.

        If ``one_to_many`` is True, allow samples to collapse into multiple
        bins if the metadata describe a one-many relationship. Supplied
        functions must allow for iteration support over the metadata key and
        must return a tuple of (path, bin) as to describe both the path in the
        hierarchy represented and the specific bin being collapsed into. The
        uniqueness of the bin is _not_ based on the path but by the name of the
        bin.

        The metadata value for the corresponding collapsed column may include
        more (or less) information about the collapsed data. For example, if
        collapsing "FOO", and there are samples that span three associations A,
        B, and C, such that sample 1 spans A and B, sample 2 spans B and C and
        sample 3 spans A and C, the resulting table will contain three
        collapsed samples:

        - A, containing original sample 1 and 3
        - B, containing original sample 1 and 2
        - C, containing original sample 2 and 3

        If a sample maps to the same bin multiple times, it will be
        counted multiple times.

        There are two supported modes for handling one-to-many relationships
        via ``one_to_many_mode``: ``add`` and ``divide``. ``add`` will add the
        sample counts to each bin that the sample maps to, which may increase
        the total number of counts in the output table. ``divide`` will divide
        a sample's counts by the number of metadata that the sample has before
        adding the counts to each bin. This will not increase the total number
        of counts in the output table.

        If ``one_to_many_md_key`` is specified, that becomes the metadata
        key that describes the collapsed path. If a value is not specified,
        then it defaults to 'Path'.

        If ``strict`` is specified, then all metadata pathways operated on
        must be indexable by ``metadata_f``.

        ``one_to_many`` and ``norm`` are not supported together.

        ``one_to_many`` and ``reduce_f`` are not supported together.

        ``one_to_many`` and ``min_group_size`` are not supported together.

        A final note on space consumption. At present, the ``one_to_many``
        functionality requires a temporary dense matrix representation. This
        was done so as it initially seems like true support requires rapid
        ``__setitem__`` functionality on the ``SparseObj`` and at the time of
        implementation, ``CSMat`` was O(N) to the number of nonzero elements.
        This is a work around until either a better ``__setitem__``
        implementation is in play on ``CSMat`` or a hybrid solution that allows
        for multiple ``SparseObj`` types is used.
        """
        if constructor is None:
            constructor = self.__class__

        collapsed_data = []
        collapsed_sample_ids = []

        if include_collapsed_metadata:
            collapsed_sample_md = []
        else:
            collapsed_sample_md = None

        if one_to_many_mode not in ['add', 'divide']:
            raise ValueError("Unrecognized one-to-many mode '%s'. Must be "
                             "either 'add' or 'divide'." % one_to_many_mode)

        if one_to_many:
            if norm:
                raise AttributeError(
                    "norm and one_to_many are not supported together")

            # determine the collapsed pathway
            # we drop all other associated metadata
            new_s_md = {}
            s_md_count = {}
            for id_, md in zip(self.sample_ids, self.sample_metadata):
                md_iter = metadata_f(md)
                num_md = 0
                while True:
                    try:
                        pathway, bin = md_iter.next()
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

                    new_s_md[bin] = pathway
                    num_md += 1

                s_md_count[id_] = num_md

            n_s = len(new_s_md)
            s_idx = dict([(bin_, i) for i, bin_ in
                         enumerate(sorted(new_s_md))])

            # We need to store floats, not ints, as things won't always divide
            # evenly.
            if one_to_many_mode == 'divide':
                dtype = float
            else:
                dtype = self._dtype

            # allocate new data. using a dense representation allows for a
            # workaround on CSMat.__setitem__ O(N) lookup. Assuming the number
            # of collapsed samples is reasonable, then this doesn't suck too
            # much.
            new_data = zeros((len(self.observation_ids), n_s), dtype=dtype)

            # for each sample
            # for each bin in the metadata
            # for each value associated with the sample
            for s_v, s_id, s_md in self.iter_samples():
                md_iter = metadata_f(s_md)
                while True:
                    try:
                        pathway, bin = md_iter.next()
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
                        new_data[:, s_idx[bin]] += s_v
                    elif one_to_many_mode == 'divide':
                        new_data[:, s_idx[bin]] += s_v / s_md_count[s_id]
                    else:
                        # Should never get here.
                        raise ValueError("Unrecognized one-to-many mode '%s'. "
                                         "Must be either 'add' or 'divide'." %
                                         one_to_many_mode)

            if include_collapsed_metadata:
                # reassociate pathway information
                for k, i in sorted(s_idx.items(), key=itemgetter(1)):
                    collapsed_sample_md.append(
                        {one_to_many_md_key: new_s_md[k]})

            # get the new sample IDs
            collapsed_sample_ids = [k for k, i in sorted(s_idx.items(),
                                                         key=itemgetter(1))]

            # convert back to self type
            data = self._conv_to_self_type(new_data)
        else:
            for bin, table in self.bin_samples_by_metadata(metadata_f):
                if len(table.sample_ids) < min_group_size:
                    continue

                redux_data = table.reduce(reduce_f, 'observation')
                if norm:
                    redux_data /= len(table.sample_ids)

                collapsed_data.append(self._conv_to_self_type(redux_data))
                collapsed_sample_ids.append(bin)

                if include_collapsed_metadata:
                    # retain metadata but store by original sample id
                    tmp_md = {}
                    for id_, md in zip(table.sample_ids,
                                       table.sample_metadata):
                        tmp_md[id_] = md
                    collapsed_sample_md.append(tmp_md)

            data = self._conv_to_self_type(collapsed_data, transpose=True)

        # if the table is empty
        if 0 in data.shape:
            raise TableException("Collapsed table is empty!")

        return table_factory(data, collapsed_sample_ids,
                             self.observation_ids[:], collapsed_sample_md,
                             self.observation_metadata, self.table_id,
                             constructor=constructor)

    def collapse_observations_by_metadata(self, metadata_f, reduce_f=add,
                                          norm=True, min_group_size=2,
                                          include_collapsed_metadata=True,
                                          constructor=None, one_to_many=False,
                                          one_to_many_mode='add',
                                          one_to_many_md_key='Path',
                                          strict=False):
        """Collapse observations in a table by observation metadata

        Bin observations by metadata then collapse each bin into a single
        observation.

        If ``include_collapsed_metadata`` is True, metadata for the collapsed
        observations are retained and can be referred to by the
        ``ObservationId`` from each observation within the bin.

        ``constructor``: the type of collapsed table to create, e.g.
        SparseTaxonTable. If None, the collapsed table will be the same type as
        this table.

        The remainder is only relevant to setting ``one_to_many`` to True.

        If ``one_to_many`` is True, allow observations to fall into multiple
        bins if the metadata describe a one-many relationship. Supplied
        functions must allow for iteration support over the metadata key and
        must return a tuple of (path, bin) as to describe both the path in the
        hierarchy represented and the specific bin being collapsed into. The
        uniqueness of the bin is _not_ based on the path but by the name of the
        bin.

        The metadata value for the corresponding collapsed row may include more
        (or less) information about the collapsed data. For example, if
        collapsing "KEGG Pathways", and there are observations that span three
        pathways A, B, and C, such that observation 1 spans A and B,
        observation 2 spans B and C and observation 3 spans A and C, the
        resulting table will contain three collapsed observations:

        - A, containing original observation 1 and 3
        - B, containing original observation 1 and 2
        - C, containing original observation 2 and 3

        If an observation maps to the same bin multiple times, it will be
        counted multiple times.

        There are two supported modes for handling one-to-many relationships
        via ``one_to_many_mode``: ``add`` and ``divide``. ``add`` will add the
        observation counts to each bin that the observation maps to, which may
        increase the total number of counts in the output table. ``divide``
        will divide an observation's counts by the number of metadata that the
        observation has before adding the counts to each bin. This will not
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
        functionality requires a temporary dense matrix representation. This
        was done so as it initially seems like true support requires rapid
        ``__setitem__`` functionality on the ``SparseObj`` and at the time of
        implementation, ``CSMat`` was O(N) to the number of nonzero elements.
        This is a work around until either a better ``__setitem__``
        implementation is in play on ``CSMat`` or a hybrid solution that allows
        for multiple ``SparseObj`` types is used.
        """
        if constructor is None:
            constructor = self.__class__

        collapsed_data = []
        collapsed_obs_ids = []

        if include_collapsed_metadata:
            collapsed_obs_md = []
        else:
            collapsed_obs_md = None

        if one_to_many_mode not in ['add', 'divide']:
            raise ValueError("Unrecognized one-to-many mode '%s'. Must be "
                             "either 'add' or 'divide'." % one_to_many_mode)

        if one_to_many:
            if norm:
                raise AttributeError(
                    "norm and one_to_many are not supported together")

            # determine the collapsed pathway
            # we drop all other associated metadata
            new_obs_md = {}
            obs_md_count = {}
            for id_, md in zip(self.observation_ids,
                               self.observation_metadata):
                md_iter = metadata_f(md)
                num_md = 0
                while True:
                    try:
                        pathway, bin = md_iter.next()
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

                    # keyed by last field in hierarchy
                    new_obs_md[bin] = pathway
                    num_md += 1

                obs_md_count[id_] = num_md

            n_obs = len(new_obs_md)
            obs_idx = dict([(bin_, i)
                           for i, bin_ in enumerate(sorted(new_obs_md))])

            # We need to store floats, not ints, as things won't always divide
            # evenly.
            if one_to_many_mode == 'divide':
                dtype = float
            else:
                dtype = self._dtype

            # allocate new data. using a dense representation allows for a
            # workaround on CSMat.__setitem__ O(N) lookup. Assuming the number
            # of collapsed observations is reasonable, then this doesn't suck
            # too much.
            new_data = zeros((n_obs, len(self.sample_ids)), dtype=dtype)

            # for each observation
            # for each bin in the metadata
            # for each value associated with the observation
            for obs_v, obs_id, obs_md in self.iter_observations():
                md_iter = metadata_f(obs_md)
                while True:
                    try:
                        pathway, bin = md_iter.next()
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
                        new_data[obs_idx[bin], :] += obs_v
                    elif one_to_many_mode == 'divide':
                        new_data[obs_idx[bin], :] += \
                            obs_v / obs_md_count[obs_id]
                    else:
                        # Should never get here.
                        raise ValueError("Unrecognized one-to-many mode '%s'. "
                                         "Must be either 'add' or 'divide'." %
                                         one_to_many_mode)

            if include_collapsed_metadata:
                # associate the pathways back
                for k, i in sorted(obs_idx.items(), key=itemgetter(1)):
                    collapsed_obs_md.append(
                        {one_to_many_md_key: new_obs_md[k]})

            # get the new observation IDs
            collapsed_obs_ids = [k for k, i in sorted(obs_idx.items(),
                                                      key=itemgetter(1))]

            # convert back to self type
            data = self._conv_to_self_type(new_data)
        else:
            for bin, table in self.bin_observations_by_metadata(metadata_f):
                if len(table.observation_ids) < min_group_size:
                    continue

                redux_data = table.reduce(reduce_f, 'sample')
                if norm:
                    redux_data /= len(table.observation_ids)

                collapsed_data.append(self._conv_to_self_type(redux_data))
                collapsed_obs_ids.append(bin)

                if include_collapsed_metadata:
                    # retain metadata but store by original observation id
                    tmp_md = {}
                    for id_, md in zip(table.observation_ids,
                                       table.observation_metadata):
                        tmp_md[id_] = md
                    collapsed_obs_md.append(tmp_md)

            data = self._conv_to_self_type(collapsed_data)

        # if the table is empty
        if 0 in data.shape:
            raise TableException("Collapsed table is empty!")

        return table_factory(data, self.sample_ids[:], collapsed_obs_ids,
                             self.sample_metadata, collapsed_obs_md,
                             self.table_id,
                             constructor=constructor)

    def transform_samples(self, f):
        """Iterate over samples, applying a function ``f`` to each value

        ``f`` must take three values: a sample value (int or float), a sample
        id, and a sample metadata entry, and return a single value (int or
        float) that replaces the provided sample value
        """
        new_m = []

        for s_v, s_id, s_md in self.iter_samples():
            new_m.append(self._conv_to_self_type(f(s_v, s_id, s_md)))

        return self.__class__(self._conv_to_self_type(new_m, transpose=True),
                              self.sample_ids[:], self.observation_ids[
                                  :], self.sample_metadata,
                              self.observation_metadata, self.table_id)

    def transform_observations(self, f):
        """Iterate over observations, applying a function ``f`` to each value

        ``f`` must take three values: an observation value (int or float), an
        observation id, and an observation metadata entry, and return a single
        value (int or float) that replaces the provided observation value

        """
        new_m = []

        for obs_v, obs_id, obs_md in self.iter_observations():
            new_m.append(self._conv_to_self_type(f(obs_v, obs_id, obs_md)))

        return self.__class__(
            self._conv_to_self_type(new_m), self.sample_ids[:],
            self.observation_ids[:], self.sample_metadata,
            self.observation_metadata, self.table_id)

    def norm_observation_by_sample(self):
        """Return new table with vals as relative abundances within each sample
        """
        def f(samp_v, samp_id, samp_md):
            return samp_v / float(samp_v.sum())
        return self.transform_samples(f)

    def norm_sample_by_observation(self):
        """Return new table with vals as relative abundances within each obs
        """
        def f(obs_v, obs_id, obs_md):
            return obs_v / float(obs_v.sum())
        # f = lambda x: x / float(x.sum())
        return self.transform_observations(f)

    def norm_observation_by_metadata(self, obs_metadata_id):
        """Return new table with vals divided by obs_metadata_id
        """
        def f(obs_v, obs_id, obs_md):
            return obs_v / obs_md[obs_metadata_id]
        return self.transform_observations(f)

    def nonzero(self):
        """Returns locations of nonzero elements within the data matrix

        The values returned are ``(observation_id, sample_id)``
        """
        # this is naively implemented. If performance is a concern, private
        # methods can be written to hit against the underlying types directly
        for o_idx, samp_vals in enumerate(self.iter_observation_data()):
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
            dtype = self._data.dtype
            op = lambda x: x.sum()

        if axis is 'sample':
            # can use np.bincount for CSMat or ScipySparse
            result = zeros(len(self.sample_ids), dtype=dtype)
            for idx, vals in enumerate(self.iter_sample_data()):
                result[idx] = op(vals)
        elif axis is 'observation':
            # can use np.bincount for CSMat or ScipySparse
            result = zeros(len(self.observation_ids), dtype=dtype)
            for idx, vals in enumerate(self.iter_observation_data()):
                result[idx] = op(vals)
        else:
            result = zeros(1, dtype=dtype)
            for vals in self.iter_sample_data():
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
            if self.sample_metadata is None or not self.sample_exists(id_):
                self_md = None
            else:
                self_md = self.sample_metadata[self_samp_idx[id_]]

            # if we have sample metadata, grab it
            if other.sample_metadata is None or not other.sample_exists(id_):
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
               not self.observation_exists(id_):
                self_md = None
            else:
                self_md = self.observation_metadata[self_obs_idx[id_]]

            # if we have observation metadata, grab it
            if other.observation_metadata is None or \
                    not other.observation_exists(id_):
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
            if other.observation_exists(obs_id):
                other_vec = other.observation_data(obs_id)
            else:
                other_vec = None

            # see if the observation exists in self, if so, pull it out.
            # if not, set to the placeholder missing
            if self.observation_exists(obs_id):
                self_vec = self.observation_data(obs_id)
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

        return self.__class__(self._conv_to_self_type(vals), sample_ids[:],
                              obs_ids[:], sample_md, obs_md)

    def format_hdf5(self, h5grp, generated_by):
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

        See Also
        --------
        Table.format_hdf5

        Examples
        --------
        ### is it okay to actually create files in doctest?

        """
        def axis_dump(grp, ids, md, order):
            """Store for an axis"""
            self._data.convert(order)

            len_ids = len(ids)
            len_indptr = len(self._data._matrix.indptr)
            len_data = self._data.size

            grp.create_dataset('data', shape=(len_data,),
                               dtype=np.float64,
                               data=self._data._matrix.data)
            grp.create_dataset('indices', shape=(len_data,),
                               dtype=np.int32,
                               data=self._data._matrix.indices)
            grp.create_dataset('indptr', shape=(len_indptr,),
                               dtype=np.int32,
                               data=self._data._matrix.indptr)

            # if we store IDs in the table as numpy arrays then this store
            # is cleaner, as is the parse
            grp.create_dataset('ids', shape=(len_ids,),
                               dtype=H5PY_VLEN_STR,
                               data=[str(i) for i in ids])

            if md is not None:
                md_str = empty(shape=(), dtype=object)
                md_str[()] = dumps(md)
                grp.create_dataset('metadata', shape=(1,),
                                   dtype=H5PY_VLEN_STR,
                                   data=md_str)

        h5grp.attrs['id'] = self.table_id if self.table_id else "No Table ID"
        h5grp.attrs['type'] = self.type
        h5grp.attrs['format-url'] = "http://biom-format.org"
        h5grp.attrs['format-version'] = (2, 0)
        h5grp.attrs['generated-by'] = generated_by
        h5grp.attrs['creation-date'] = datetime.now().isoformat()
        h5grp.attrs['shape'] = self._data.shape
        h5grp.attrs['nnz'] = self._data.size

        axis_dump(h5grp.create_group('observation'), self.observation_ids,
                  self.observation_metadata, 'csr')
        axis_dump(h5grp.create_group('sample'), self.sample_ids,
                  self.sample_metadata, 'csc')

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
            num_rows, num_cols = self._data.shape
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
        for obs_index, obs in enumerate(self.iter_observations()):
            biom_format_obj["rows"].append(
                {"id": "%s" % obs[1], "metadata": obs[2]})
            # If the matrix is dense, simply convert the numpy array to a list
            # of data values. If the matrix is sparse, we need to store the
            # data in sparse format, as it is given to us in a numpy array in
            # dense format (i.e. includes zeroes) by iter_observations().
            dense_values = list(obs[0])
            sparse_values = []
            for col_index, val in enumerate(dense_values):
                if float(val) != 0.0:
                    sparse_values.append([obs_index, col_index,
                                          dtype(val)])
            biom_format_obj["data"].extend(sparse_values)

        # Fill in details about the columns in the table.
        biom_format_obj["columns"] = []
        for samp in self.iter_samples():
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
            num_rows, num_cols = self._data.shape
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
        for obs_index, obs in enumerate(self.iter_observations()):
            # i'm crying on the inside
            if obs_index != max_row_idx:
                rows.append('{"id": "%s", "metadata": %s},' % (obs[1],
                                                               dumps(obs[2])))
            else:
                rows.append('{"id": "%s", "metadata": %s}],' % (obs[1],
                                                                dumps(obs[2])))

            # if we are not on the last row
            # if obs_index != max_row_idx:
            #    data.append("[%s]," % ','.join(map(repr, obs[0])))
            # else:
            #    data.append("[%s]]," % ','.join(map(repr, obs[0])))

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
        for samp_index, samp in enumerate(self.iter_samples()):
            if samp_index != max_col_idx:
                columns.append('{"id": "%s", "metadata": %s},' % (samp[1],
                                                                  dumps(
                                                                  samp[2])))
            else:
                columns.append('{"id": "%s", "metadata": %s}]' % (samp[1],
                                                                  dumps(
                                                                  samp[2])))

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
        for (foo, col_idx), val in row.iteritems():
            mat[row_idx, col_idx] = val

    return mat


def table_factory(data, sample_ids, observation_ids, sample_metadata=None,
                  observation_metadata=None, table_id=None, **kwargs):
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
                          sample_ids,
                          observation_ids,
                          sample_md,
                          observation_md)
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
        data = nparray_to_sparseobj(data, dtype)

    # if we have a list of things
    elif isinstance(data, list):
        if not data:
            raise TableException("No data was supplied. Cannot create "
                                 "an empty table.")

        elif isinstance(data[0], ndarray):
            data = list_nparray_to_sparseobj(data, dtype)

        elif isinstance(data[0], dict):
            data = list_dict_to_sparseobj(data, dtype)

        elif isinstance(data[0], list):
            data = list_list_to_sparseobj(data, dtype, shape=shape)

        else:
            raise TableException("Unknown nested list type")

    # if we have a dict representation
    elif isinstance(data, dict) and not isinstance(data, SparseObj):
        data = dict_to_sparseobj(data, dtype)

    elif isinstance(data, tuple) and isinstance(data[0], ndarray):
        # give it a go...
        # there isn't a CSMat equivilent
        from biom.backends.scipysparse import coo_arrays_to_scipy
        data = coo_arrays_to_scipy(data)

    elif isinstance(data, SparseObj):
        pass

    else:
        raise TableException("Cannot handle data!")

    return Table(data, sample_ids, observation_ids,
                 sample_metadata=sample_metadata,
                 observation_metadata=observation_metadata,
                 table_id=table_id, **kwargs)


def get_zerod_matrix(mat, dtype=float):
    """Returns a zerod matrix"""
    if isinstance(mat, ndarray):
        return zeros(mat.shape, dtype=float)
    elif isinstance(mat, SparseObj):
        return SparseObj(*mat.shape, dtype=float)
    else:
        raise TableException("Unknown mat type")
