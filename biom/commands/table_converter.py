#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division
from pyqi.core.command import (Command, CommandIn, CommandOut,
                               ParameterCollection)
from pyqi.core.exception import CommandError
from biom.table import Table
from biom.parse import MetadataMap

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Daniel McDonald",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout",
               "Jose Antonio Navas Molina", "Jorge Cañardo Alastuey"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


class TableConverter(Command):
    TableTypes = ["OTU table",
                  "Pathway table",
                  "Function table",
                  "Ortholog table",
                  "Gene table",
                  "Metabolite table",
                  "Taxon table"]

    ObservationMetadataTypes = {
        'sc_separated': lambda x: [e.strip() for e in x.split(';')],
        'naive': lambda x: x
    }

    ObservationMetadataFormatters = {
        'sc_separated': lambda x: '; '.join(x),
        'naive': lambda x: x
    }

    ObservationMetadataTypes['taxonomy'] = \
        ObservationMetadataTypes['sc_separated']

    BriefDescription = "Convert to/from the BIOM table format"
    LongDescription = ("Convert between BIOM and 'classic' (tab-delimited) "
                       "table formats. Detailed usage examples can be found "
                       "here: http://biom-format.org/documentation/biom_conver"
                       "sion.html")

    CommandIns = ParameterCollection([
        # This is not an ideal usage of the pyqi framework because we are
        # expecting a file-like object here, and a lot of the parameters deal
        # with I/O-ish things, like converting between file formats. Even
        # though no I/O is forced here, it would be better to have rich objects
        # as input and output, instead of lines of data. However, this will
        # likely require a refactoring/redesign of our interface for table
        # conversions because the primary input here can be either a BIOM table
        # or a classic table. One possible solution is to split out different
        # types of conversions into their own (smaller and simpler) commands,
        # which would allow us to avoid some of this I/O-ish stuff.
        CommandIn(Name='table', DataType=Table,
                  Description='the input table (file-like object), either in '
                  'BIOM or classic format', Required=True),
        CommandIn(Name='to_json', DataType=bool,
                  Description='Output as a JSON table', Default=False),
        CommandIn(Name='to_hdf5', DataType=bool,
                  Description='Output as a HDF5 table', Default=False),
        CommandIn(Name='to_tsv', DataType=bool,
                  Description='Output as a TSV table', Default=False),
        CommandIn(Name='collapsed_samples', DataType=bool,
                  Description='If to_hdf5 and the original table is a '
                              'collapsed by samples biom table, this will '
                              'update the sample metadata of the table to '
                              'the supported HDF5 collapsed format'),
        CommandIn(Name='collapsed_observations', DataType=bool,
                  Description='If to_hdf5 and the original table is a '
                              'collapsed by observations biom table, this will'
                              ' update the observation metadata of the table '
                              'to the supported HDF5 collapsed format'),
        CommandIn(Name='sample_metadata', DataType=MetadataMap,
                  Description='the sample metadata map (will add sample '
                  'metadata to the BIOM table, if provided). Only applies '
                  'when converting from classic table file to BIOM table '
                  'file'),
        CommandIn(Name='observation_metadata', DataType=MetadataMap,
                  Description='the observation metadata map (will add '
                  'observation metadata to the BIOM table, if provided). Only '
                  'applies when converting from classic table file to BIOM '
                  'table file'),
        CommandIn(Name='header_key', DataType=str,
                  Description='pull this key from observation metadata within '
                  'a BIOM table file when creating a classic table file',
                  DefaultDescription='no observation metadata will be '
                  'included'),
        CommandIn(Name='output_metadata_id', DataType=str,
                  Description='the name to be given to the observation '
                  'metadata column when creating a classic table from a BIOM-'
                  'formatted table', DefaultDescription='same name as in the '
                  'BIOM-formatted table'),
        CommandIn(Name='table_type', DataType=str,
                  Description='the type of the table, must be one of: %s' %
                  ', '.join(TableTypes), Required=False),
        CommandIn(Name='process_obs_metadata', DataType=str,
                  Description='process metadata associated with observations '
                  'when converting from a classic table. Must be one of: %s' %
                  ', '.join(ObservationMetadataTypes), Default=None),
        CommandIn(Name='tsv_metadata_formatter', DataType=str,
                  Description='Method for formatting the observation '
                  'metadata, must be one of: %s' %
                  ', '.join(ObservationMetadataFormatters),
                  Default='sc_separated')
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name='table', DataType=tuple,
                   Description='The resulting table and format')
    ])

    def run(self, **kwargs):
        table = kwargs['table']
        sample_metadata = kwargs['sample_metadata']
        observation_metadata = kwargs['observation_metadata']
        header_key = kwargs['header_key']
        output_metadata_id = kwargs['output_metadata_id']
        process_obs_metadata = kwargs['process_obs_metadata']
        obs_md_fmt = kwargs['tsv_metadata_formatter']
        table_type = kwargs['table_type']
        to_tsv = kwargs['to_tsv']
        to_hdf5 = kwargs['to_hdf5']
        to_json = kwargs['to_json']
        collapsed_observations = kwargs['collapsed_observations']
        collapsed_samples = kwargs['collapsed_samples']

        if sum([to_tsv, to_hdf5, to_json]) == 0:
            raise CommandError("Must specify an output format")
        elif sum([to_tsv, to_hdf5, to_json]) > 1:
            raise CommandError("Can only specify a single output format")

        # if we don't have a table type, then one is required to be specified
        if table.type in [None, "None"]:
            if table_type is None:
                raise CommandError("Must specify --table-type!")
            else:
                if table_type not in self.TableTypes:
                    raise CommandError("Unknown table type: %s" % table_type)

                table.type = table_type

        if obs_md_fmt not in self.ObservationMetadataFormatters:
            raise CommandError("Unknown tsv_metadata_formatter: %s" %
                               obs_md_fmt)
        else:
            obs_md_fmt_f = self.ObservationMetadataFormatters[obs_md_fmt]

        if sample_metadata is not None:
            table.add_metadata(sample_metadata)

        # if the user does not specify a name for the output metadata column,
        # set it to the same as the header key
        output_metadata_id = output_metadata_id or header_key

        if process_obs_metadata is not None and not to_tsv:
            if process_obs_metadata not in self.ObservationMetadataTypes:
                raise CommandError(
                    "Unknown observation metadata processing method, must be "
                    "one of: %s" %
                    ', '.join(self.ObservationMetadataTypes.keys()))

            if table.metadata(axis='observation') is None:
                raise CommandError("Observation metadata processing requested "
                                   "but it doesn't appear that there is any "
                                   "metadata to operate on!")

            # and if this came in as TSV, then we expect only a single type of
            # metadata
            md_key = table.metadata(axis='observation')[0].keys()[0]

            process_f = self.ObservationMetadataTypes[process_obs_metadata]
            it = zip(table.ids(axis='observation'),
                     table.metadata(axis='observation'))
            new_md = {id_: {md_key: process_f(md[md_key])} for id_, md in it}

            if observation_metadata:
                for k, v in observation_metadata.items():
                    new_md[k].update(v)
            table.add_metadata(new_md, 'observation')

        if to_tsv:
            result = table.to_tsv(header_key=header_key,
                                  header_value=output_metadata_id,
                                  metadata_formatter=obs_md_fmt_f)
            fmt = 'tsv'
        elif to_json:
            result = table
            fmt = 'json'
        elif to_hdf5:
            result = table
            if collapsed_observations:
                metadata = [{'collapsed_ids': md.keys()}
                            for md in result.metadata(axis='observation')]
                result._observation_metadata = metadata
            if collapsed_samples:
                metadata = [{'collapsed_ids': md.keys()}
                            for md in result.metadata()]
                result._sample_metadata = metadata
            if collapsed_observations or collapsed_samples:
                # We have changed the metadata, it is safer to make sure that
                # it is correct
                result._cast_metadata()
            fmt = 'hdf5'

        return {'table': (result, fmt)}

CommandConstructor = TableConverter
