# -----------------------------------------------------------------------------
# Copyright (c) 2011-2015, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division

from biom.parse import MetadataMap


def _split_on_semicolons(x):
    return [e.strip() for e in x.split(';')]


def _split_on_semicolons_and_pipes(x):
    return [[e.strip() for e in y.split(';')] for y in x.split('|')]


def _int(x):
    try:
        return int(x)
    except ValueError:
        return x


def _float(x):
    try:
        return float(x)
    except ValueError:
        return x


def add_metadata(table, sample_metadata=None, observation_metadata=None,
                 sc_separated=None, sc_pipe_separated=None, int_fields=None,
                 float_fields=None, sample_header=None,
                 observation_header=None):

    if sample_metadata is None and observation_metadata is None:
        raise ValueError('Must specify sample_metadata and/or '
                         'observation_metadata.')

    # define metadata processing functions, if any
    process_fns = {}
    if sc_separated is not None:
        process_fns.update(dict.fromkeys(sc_separated,
                                         _split_on_semicolons))

    if sc_pipe_separated is not None:
        process_fns.update(dict.fromkeys(sc_pipe_separated,
                           _split_on_semicolons_and_pipes))

    if int_fields is not None:
        process_fns.update(dict.fromkeys(int_fields, _int))

    if float_fields is not None:
        process_fns.update(dict.fromkeys(float_fields, _float))

    # parse mapping files
    if sample_metadata is not None:
        sample_metadata = MetadataMap.from_file(sample_metadata,
                                                process_fns=process_fns,
                                                header=sample_header)

    if observation_metadata is not None:
        observation_metadata = MetadataMap.from_file(
            observation_metadata,
            process_fns=process_fns,
            header=observation_header)

    # NAUGHTY: this is modifying the input table IN PLACE!!! And then
    # RETURNING IT! MetadataAdder is angry!

    # add metadata as necessary
    if sample_metadata:
        table.add_metadata(sample_metadata, axis='sample')

    if observation_metadata:
        table.add_metadata(observation_metadata, axis='observation')

    return table
