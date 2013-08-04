#!/usr/bin/env python

from __future__ import division
from pyqi.core.command import Command, Parameter
from biom.parse import MetadataMap
from biom.table import Table

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

class AddMetadata(Command):
    BriefDescription = "Add metadata to table"
    LongDescription = ("Add sample and/or observation metadata to "
                       "BIOM-formatted files. Detailed usage examples can be "
                       "found here: http://biom-format.org/documentation/adding_metadata.html")

    def run(self, **kwargs):
        results = {}

        sc_separated = kwargs['sc-separated']
        int_fields = kwargs['int-fields']
        float_fields = kwargs['float-fields']
        sc_pipe_separated = kwargs['sc-pipe-separated']
        sample_mapping = kwargs['sample-mapping']
        obs_mapping = kwargs['observation-mapping']

        ## define metadata processing functions, if any
        process_fns = {}
        if sc_separated != None:
            process_fns.update({}.fromkeys(sc_separated,
                               self._split_on_semicolons))

        if int_fields != None:
            process_fns.update({}.fromkeys(int_fields, self._int))

        if float_fields != None:
            process_fns.update({}.fromkeys(float_fields, self._float))

        if sc_pipe_separated != None:
            process_fns.update({}.fromkeys(sc_pipe_separated,
                    self._split_on_semicolons_and_pipes))

        ## parse mapping files
        if sample_mapping != None:
            sample_mapping = MetadataMap.fromFile(sample_mapping,
                                                  process_fns=process_fns,
                                                  header=sample_header)

        if obs_mapping != None:
            obs_mapping = MetadataMap.fromFile(obs_mapping,
                                               process_fns=process_fns,
                                               header=observation_header)

        if sample_mapping == None and obs_mapping == None:
            raise CommandError('Must specify sample_mapping and/or '
                               'obs_mapping.')

        table = kwargs['biom-table']

        ## add metadata as necessary
        if sample_mapping:
            table.addSampleMetadata(sample_mapping)

        if obs_mapping:
            table.addObservationMetadata(obs_mapping)

        # NAUGHTY: this is modifying the input table IN PLACE!!! And then
        # RETURNING IT!
        results['biom-table'] = table

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

    def _get_parameters(self):
        return [Parameter(Name='biom-table', Required=True, Type=Table,
                          Help='the input BIOM table'),
                Parameter(Name='sample-mapping', Required=False,
                          Type=MetadataMap, Help='the sample metadata map '
                          '(will add sample metadata to the input BIOM table, '
                          'if provided)', Default=None),
                Parameter(Name='observation-mapping', Required=False,
                          Type=MetadataMap, Help='the observation metadata '
                          'map (will add observation metadata to the input '
                          'BIOM table, if provided)', Default=None),
                Parameter(Name='sc-separated', Required=False,
                          Type=list, Help='list of the metadata fields '
                          'to split on semicolons. This is useful for '
                          'hierarchical data such as taxonomy or functional '
                          'categories', Default=None),
                Parameter(Name='sc-pipe-separated', Required=False,
                          Type=list, Help='list of the metadata fields to '
                          'split on semicolons and pipes ("|"). This is '
                          'useful for hierarchical data such as functional '
                          'categories with one-to-many mappings (e.g. '
                          'x;y;z|x;y;w)', Default=None),
                Parameter(Name='int-fields', Required=False,
                          Type=list, Help='list of the metadata fields to '
                          'cast to integers. This is useful for integer data '
                          'such as "DaysSinceStart"', Default=None),
                Parameter(Name='float-fields', Required=False,
                          Type=list, Help='list of the metadata fields to '
                          'cast to floating point numbers. This is useful for '
                          'real number data such as "pH"', Default=None),
                Parameter(Name='observation-header', Required=False,
                          Type=list, Help='list of the observation metadata '
                          'field names. This is useful if a header line is '
                          'not provided with the metadata, if you want to '
                          'rename the fields, or if you want to include only '
                          'the first n fields where n is the number of '
                          'entries provided here', Default=None,
                          DefaultDescription='use header from observation-'
                          'mapping'),
                Parameter(Name='sample-header', Required=False,
                          Type=list, Help='list of the sample metadata field '
                          'names. This is useful if a header line is not '
                          'provided with the metadata, if you want to rename '
                          'the fields, or if you want to include only the '
                          'first n fields where n is the number of entries '
                          'provided here', Default=None,
                          DefaultDescription='use header from sample-mapping')]

CommandConstructor = AddMetadata
