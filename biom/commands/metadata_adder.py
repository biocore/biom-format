#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from pyqi.core.command import (Command, CommandIn, CommandOut, 
        ParameterCollection)
from pyqi.core.exception import CommandError
from biom.parse import MetadataMap
from biom.table import Table

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Morgan Langille", "Jai Ram Rideout",
               "Daniel McDonald"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class MetadataAdder(Command):
    BriefDescription = "Add metadata to a BIOM table"
    LongDescription = ("Add sample and/or observation metadata to "
                       "BIOM-formatted files. Detailed usage examples can be "
                       "found here: http://biom-format.org/documentation/adding_metadata.html")

    CommandIns = ParameterCollection([
        CommandIn(Name='table', DataType=Table,
                  Description='the input BIOM table', Required=True),
        # sample_metadata and observation_metadata are currently files (or
        # file-like) because of the existing metadata map / processing function
        # support. Ideally, these two parameters should be MetadataMap
        # instances.
        CommandIn(Name='sample_metadata', DataType=file,
                  Description='the sample metadata map (will add sample '
                  'metadata to the input BIOM table, if provided)'),
        CommandIn(Name='observation_metadata', DataType=file,
                  Description='the observation metadata map (will add '
                  'observation metadata to the input BIOM table, if '
                  'provided)'),
        CommandIn(Name='sc_separated', DataType=list,
                  Description='list of the metadata fields to split on '
                  'semicolons. This is useful for hierarchical data such as '
                  'taxonomy or functional categories'),
        CommandIn(Name='sc_pipe_separated', DataType=list,
                  Description='list of the metadata fields to split on '
                  'semicolons and pipes ("|"). This is useful for '
                  'hierarchical data such as functional categories with '
                  'one-to-many mappings (e.g. x;y;z|x;y;w)'),
        CommandIn(Name='int_fields', DataType=list,
                  Description='list of the metadata fields to cast to '
                  'integers. This is useful for integer data such as '
                  '"DaysSinceStart"'),
        CommandIn(Name='float_fields', DataType=list,
                  Description='list of the metadata fields to cast to '
                  'floating point numbers. This is useful for real number '
                  'data such as "pH"'),
        CommandIn(Name='sample_header', DataType=list,
                  Description='list of the sample metadata field names. This '
                  'is useful if a header line is not provided with the '
                  'metadata, if you want to rename the fields, or if you want '
                  'to include only the first n fields where n is the number '
                  'of entries provided here',
                  DefaultDescription='use header from sample metadata map'),
        CommandIn(Name='observation_header', DataType=list,
                  Description='list of the observation metadata field names. '
                  'This is useful if a header line is not provided with the '
                  'metadata, if you want to rename the fields, or if you want '
                  'to include only the first n fields where n is the number '
                  'of entries provided here',
                  DefaultDescription='use header from observation metadata '
                  'map')
    ])

    CommandOuts = ParameterCollection([
            CommandOut(Name='table', DataType=Table,
                       Description='Table with added metadata')
            ])

    def run(self, **kwargs):
        table = kwargs['table']
        sample_metadata = kwargs['sample_metadata']
        observation_metadata = kwargs['observation_metadata']
        sc_separated = kwargs['sc_separated']
        sc_pipe_separated = kwargs['sc_pipe_separated']
        int_fields = kwargs['int_fields']
        float_fields = kwargs['float_fields']
        sample_header = kwargs['sample_header']
        observation_header = kwargs['observation_header']

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
            sample_metadata = MetadataMap.fromFile(sample_metadata,
                                                   process_fns=process_fns,
                                                   header=sample_header)

        if observation_metadata is not None:
            observation_metadata = MetadataMap.fromFile(observation_metadata,
                    process_fns=process_fns, header=observation_header)

        if sample_metadata is None and observation_metadata is None:
            raise CommandError('Must specify sample_metadata and/or '
                               'observation_metadata.')

        # NAUGHTY: this is modifying the input table IN PLACE!!! And then
        # RETURNING IT! MetadataAdder is angry!

        # add metadata as necessary
        if sample_metadata:
            table.addSampleMetadata(sample_metadata)

        if observation_metadata:
            table.addObservationMetadata(observation_metadata)

        return {'table': table}

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

CommandConstructor = MetadataAdder
