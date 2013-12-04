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
from biom.table import (SparseOTUTable, DenseOTUTable, SparsePathwayTable,
                        DensePathwayTable, SparseFunctionTable,
                        DenseFunctionTable, SparseOrthologTable,
                        DenseOrthologTable, SparseGeneTable, DenseGeneTable,
                        SparseMetaboliteTable, DenseMetaboliteTable,
                        SparseTaxonTable, DenseTaxonTable, table_factory)
from biom.parse import (parse_biom_table, MetadataMap, convert_biom_to_table,
                        convert_table_to_biom, generatedby)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Daniel McDonald",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class TableConverter(Command):
    MatrixTypes = ['sparse', 'dense']

    TableTypes = {
            'otu table': [SparseOTUTable, DenseOTUTable],
            'pathway table': [SparsePathwayTable, DensePathwayTable],
            'function table': [SparseFunctionTable, DenseFunctionTable],
            'ortholog table': [SparseOrthologTable, DenseOrthologTable],
            'gene table': [SparseGeneTable, DenseGeneTable],
            'metabolite table': [SparseMetaboliteTable, DenseMetaboliteTable],
            'taxon table': [SparseTaxonTable, DenseTaxonTable]
    }

    ObservationMetadataTypes = {
            'sc_separated': lambda x: [e.strip() for e in x.split(';')],
            'naive': lambda x: x
    }
    ObservationMetadataTypes['taxonomy'] = \
            ObservationMetadataTypes['sc_separated']

    BriefDescription = "Convert to/from the BIOM table format"
    LongDescription = ("Convert between BIOM and 'classic' (tab-delimited) "
                       "table formats. Detailed usage examples can be found "
                       "here: http://biom-format.org/documentation/biom_conversion.html")

    CommandIns = ParameterCollection([
        ### This is not an ideal usage of the pyqi framework because we are
        # expecting a file-like object here, and a lot of the parameters deal
        # with I/O-ish things, like converting between file formats. Even
        # though no I/O is forced here, it would be better to have rich objects
        # as input and output, instead of lines of data. However, this will
        # likely require a refactoring/redesign of our interface for table
        # conversions because the primary input here can be either a BIOM table
        # or a classic table. One possible solution is to split out different
        # types of conversions into their own (smaller and simpler) commands,
        # which would allow us to avoid some of this I/O-ish stuff.
        CommandIn(Name='table_file', DataType=file,
                  Description='the input table (file-like object), either in '
                  'BIOM or classic format', Required=True),
        CommandIn(Name='matrix_type', DataType=str,
                  Description='the type of BIOM file to create (dense or '
                  'sparse) when a classic table is supplied',
                  Default='sparse'),
        CommandIn(Name='biom_to_classic_table', DataType=bool,
                  Description='convert BIOM table file to classic table file',
                  Default=False, DefaultDescription='convert classic table '
                  'file to BIOM table file'),
        CommandIn(Name='sparse_biom_to_dense_biom', DataType=bool,
                  Description='convert sparse BIOM table file to a dense BIOM '
                  'table file', Default=False, DefaultDescription='convert '
                  'classic table file to BIOM table file'),
        CommandIn(Name='dense_biom_to_sparse_biom', DataType=bool,
                  Description='convert dense BIOM table file to a sparse BIOM '
                  'table file', Default=False, DefaultDescription='convert '
                  'classic table file to BIOM table file'),
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
        CommandIn(Name='process_obs_metadata', DataType=str,
                  Description='process metadata associated with observations '
                  'when converting from a classic table. Must be one of: %s' %
                  ', '.join(ObservationMetadataTypes.keys()), Default='naive'),
        CommandIn(Name='table_type', DataType=str,
                  Description='the BIOM table type to get converted into. '
                  'Required when converting a classic table file to a BIOM '
                  'table file. Must be one of: %s' %
                  ', '.join(TableTypes.keys()))
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name='table_str', DataType=str,
                   Description='The resulting table')
    ])

    def run(self, **kwargs):
        table_file = kwargs['table_file']
        matrix_type = kwargs['matrix_type']
        biom_to_classic_table = kwargs['biom_to_classic_table']
        sparse_biom_to_dense_biom = kwargs['sparse_biom_to_dense_biom']
        dense_biom_to_sparse_biom = kwargs['dense_biom_to_sparse_biom']
        sample_metadata = kwargs['sample_metadata']
        observation_metadata = kwargs['observation_metadata']
        header_key = kwargs['header_key']
        output_metadata_id = kwargs['output_metadata_id']
        process_obs_metadata = kwargs['process_obs_metadata']
        table_type = kwargs['table_type']

        if sum([biom_to_classic_table, sparse_biom_to_dense_biom,
                dense_biom_to_sparse_biom]) > 1:
            raise CommandError("Converting between classic/BIOM formats and "
                               "sparse/dense representations are mutually "
                               "exclusive. You may only specify a single "
                               "operation at a time.")

        # if the user does not specify a name for the output metadata column,
        # set it to the same as the header key
        output_metadata_id = output_metadata_id or header_key

        convert_error_msg = ("Input does not look like a BIOM-formatted file. "
                             "Did you accidentally specify that a classic "
                             "table file should be created from a BIOM table "
                             "file?")
        if biom_to_classic_table:
            try:
                result = convert_biom_to_table(table_file, header_key,
                                               output_metadata_id)
            except ValueError:
                raise CommandError(convert_error_msg)
        elif sparse_biom_to_dense_biom:
            try:
                table = parse_biom_table(table_file)
            except ValueError:
                raise CommandError(convert_error_msg)

            conv_constructor = self.TableTypes[table._biom_type.lower()][1]
            conv_table = table_factory(table._data, table.SampleIds,
                            table.ObservationIds, table.SampleMetadata,
                            table.ObservationMetadata, table.TableId,
                            constructor=conv_constructor)
            result = conv_table.getBiomFormatJsonString(generatedby())
        elif dense_biom_to_sparse_biom:
            try:
                table = parse_biom_table(table_file)
            except ValueError:
                raise CommandError(convert_error_msg)

            conv_constructor = self.TableTypes[table._biom_type.lower()][0]
            conv_table = table_factory(table._data, table.SampleIds, 
                            table.ObservationIds, table.SampleMetadata, 
                            table.ObservationMetadata, table.TableId, 
                            constructor=conv_constructor)
            result = conv_table.getBiomFormatJsonString(generatedby())
        else:
            if table_type is None:
                raise CommandError("Must specify the BIOM table type: %s" %
                                   ', '.join(self.TableTypes.keys()))
            else:
                table_type = table_type.lower()

            if table_type not in self.TableTypes:
                raise CommandError("Unknown BIOM table type, must be one of: "
                                   "%s" % ', '.join(self.TableTypes.keys()))

            if matrix_type not in self.MatrixTypes:
                raise CommandError("Unknown BIOM matrix type, must be one of: "
                                   "%s" % ', '.join(self.MatrixTypes))

            if process_obs_metadata not in \
                    self.ObservationMetadataTypes.keys():
                raise CommandError("Unknown observation metadata processing "
                        "method, must be one of: %s" %
                        ', '.join(self.ObservationMetadataTypes.keys()))

            idx = 0 if matrix_type == 'sparse' else 1
            constructor = self.TableTypes[table_type][idx]


            convert_error_msg = ("Input does not look like a classic table. "
                                 "Did you forget to specify that a classic "
                                 "table file should be created from a BIOM "
                                 "table file?")
            try:
                result = convert_table_to_biom(table_file, sample_metadata,
                        observation_metadata,
                        self.ObservationMetadataTypes[process_obs_metadata],
                        constructor)
            except ValueError:
                raise CommandError(convert_error_msg)
            except IndexError:
                raise CommandError(convert_error_msg)

        return {'table_str': result}

CommandConstructor = TableConverter
