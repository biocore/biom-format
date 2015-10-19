# -----------------------------------------------------------------------------
# Copyright (c) 2011-2015, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division

from biom.parse import MetadataMap

def add_metadata(table, sample_metadata, observation_metadata,
                 sc_separated, sc_pipe_separated, int_fields, float_fields,
                 sample_header, observation_header):

    # define metadata processing functions, if any
    process_fns = {}
    if sc_separated is not None:
        process_fns.update(dict.fromkeys(sc_separated,
                                         self._split_on_semicolons))

    if sc_pipe_separated is not None:
        process_fns.update(dict.fromkeys(sc_pipe_separated,
                           self._split_on_semicolons_and_pipes))

    if int_fields is not None:
        process_fns.update(dict.fromkeys(int_fields, self._int))

    if float_fields is not None:
        process_fns.update(dict.fromkeys(float_fields, self._float))

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

    if sample_metadata is None and observation_metadata is None:
        raise CommandError('Must specify sample_metadata and/or '
                           'observation_metadata.')

    # NAUGHTY: this is modifying the input table IN PLACE!!! And then
    # RETURNING IT! MetadataAdder is angry!

    # add metadata as necessary
    if sample_metadata:
        table.add_metadata(sample_metadata, axis='sample')

    if observation_metadata:
        table.add_metadata(observation_metadata, axis='observation')

    return table

    def _split_on_semicolons(self, x):
        return [e.strip() for e in x.split(';')]

    def _split_on_semicolons_and_pipes(self, x):
        return [[e.strip() for e in y.split(';')] for y in x.split('|')]

    def _int(self, x):
        try:
            return int(x)
        except ValueError:
            return x

    def _float(self, x):
        try:
            return float(x)
        except ValueError:
            return x
