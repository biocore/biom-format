#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
import json
import dateutil.parser
from operator import and_
from pyqi.core.command import (Command, CommandIn, CommandOut, 
    ParameterCollection)

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jose Clemente", "Greg Caporaso",
               "Jai Ram Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", "Sue Huse"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__author__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class TableValidator(Command):
    BriefDescription = "Validate a BIOM-formatted file"
    LongDescription = ("Test a file for adherence to the Biological "
                       "Observation Matrix (BIOM) format specification. This "
                       "specification is defined at http://biom-format.org")

    CommandIns = ParameterCollection([
        CommandIn(Name='table_json', DataType=dict,
                  Description='the input BIOM JSON object (e.g., the output '
                  'of json.load)', Required=True),
        CommandIn(Name='format_version', DataType=str,
                  Description='the specific format version to validate '
                  'against', Required=False,
                  Default='Biological Observation Matrix 1.0.0'),
        CommandIn(Name='detailed_report', DataType=bool,
                  Description='include more details in the output report',
                  Required=False, Default=False)
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name='valid_table',
                   Description='Is the table valid?',
                   DataType=bool),
        CommandOut(Name='report_lines',
                   Description='Detailed report',
                   DataType=list)
    ])

    FormatURL = "http://biom-format.org"
    TableTypes = set(['otu table', 'pathway table', 'function table',
                      'ortholog table', 'gene table', 'metabolite table',
                      'taxon table'])
    MatrixTypes = set(['sparse', 'dense'])
    ElementTypes = {'int': int,'str': str,'float': float, 'unicode': unicode}

    def run(self, **kwargs):
        table_json = kwargs['table_json']
        # Need to make this an attribute so that we have this info during
        # validation.
        self._format_version = kwargs['format_version']
        detailed_report = kwargs['detailed_report']

        report_lines = []
        valid_table = True

        if detailed_report:
            report_lines.append("Validating BIOM table...")

        required_keys = [
                ('format', self._valid_format),
                ('format_url', self._valid_format_url),
                ('type', self._valid_type),
                ('rows', self._valid_rows),
                ('columns', self._valid_columns),
                ('shape', self._valid_shape),
                ('data', self._valid_data),
                ('matrix_type', self._valid_matrix_type),
                ('matrix_element_type', self._valid_matrix_element_type),
                ('generated_by', self._valid_generated_by),
                ('id', self._valid_nullable_id),
                ('date', self._valid_datetime)
        ]

        for key, method in required_keys:
            if key not in table_json:
                valid_table = False
                report_lines.append("Missing field: '%s'" % key)
                continue

            if detailed_report:
                report_lines.append("Validating '%s'..." % key)

            status_msg = method(table_json)

            if len(status_msg) > 0:
                valid_table = False
                report_lines.append(status_msg)

        if 'shape' in table_json:
            if detailed_report:
                report_lines.append("Validating 'shape' versus number of rows "
                                    "and columns...")

            if ('rows' in table_json and
                len(table_json['rows']) != table_json['shape'][0]):
                valid_table = False
                report_lines.append("Number of rows in 'rows' is not equal to "
                                    "'shape'")

            if ('columns' in table_json and
                len(table_json['columns']) != table_json['shape'][1]):
                valid_table = False
                report_lines.append("Number of columns in 'columns' is not "
                                    "equal to 'shape'")

        return {'valid_table': valid_table, 'report_lines': report_lines}

    def _is_int(self, x):
        """Return True if x is an int"""
        return isinstance(x, int)

    def _valid_format_url(self, table_json):
        """Check if format_url is correct"""
        if table_json['format_url'] != self.FormatURL:
            return "Invalid 'format_url'"
        else:
            return ''

    def _valid_shape(self, table_json):
        """A matrix header is (int, int) representing the size of a 2D matrix"""
        a,b = table_json['shape']

        if not (self._is_int(a) and self._is_int(b)):
            return "'shape' values do not appear to be integers"
        else:
            return ''

    def _valid_matrix_type(self, table_json):
        """Check if a valid matrix type exists"""
        if table_json['matrix_type'] not in self.MatrixTypes:
            return "Unknown 'matrix_type'"
        else:
            return ''

    def _valid_matrix_element_type(self, table_json):
        """Check if a valid element type exists"""
        if table_json['matrix_element_type'] not in self.ElementTypes:
            return "Unknown 'matrix_element_type'"
        else:
            return ''

    def _valid_datetime(self, table_json):
        """Verify datetime can be parsed

        Expects ISO 8601 datetime format (for example, 2011-12-19T19:00:00
                                          note that a 'T' separates the date 
                                          and time)
        """
        try:
            _ = dateutil.parser.parse(table_json['date'])
        except:
            return "Timestamp does not appear to be ISO 8601"
        else:
            return ''

    def _valid_sparse_data(self, table_json):
        """All index positions must be integers and values are of dtype"""
        dtype = self.ElementTypes[table_json['matrix_element_type']]
        n_rows, n_cols = table_json['shape']
        n_rows -= 1 # adjust for 0-based index
        n_cols -= 1 # adjust for 0-based index

        for idx, coord in enumerate(table_json['data']):
            try:
                x,y,val = coord
            except:
                return "Bad matrix entry idx %d: %s" % (idx,repr(coord))

            if not self._is_int(x) or not self._is_int(y):
                return "Bad x or y type at idx %d: %s" % (idx,repr(coord))

            if not isinstance(val, dtype):
                return "Bad value at idx %d: %s" % (idx,repr(coord))

            if x < 0 or x > n_rows:
                return "x out of bounds at idx %d: %s" % (idx,repr(coord))

            if y < 0 or y > n_cols:
                return "y out of bounds at idx %d: %s" % (idx,repr(coord))

        return ''

    def _valid_dense_data(self, table_json):
        """All elements must be of dtype and correspond to shape"""
        dtype = self.ElementTypes[table_json['matrix_element_type']]
        n_rows, n_cols = table_json['shape']

        for row in table_json['data']:
            if len(row) != n_cols:
                return "Incorrect number of cols: %s" % repr(row)

            if not reduce(and_, [isinstance(v, dtype) for v in row]):
                return "Bad datatype in row: %s" % repr(row)

        if len(table_json['data']) != n_rows:
            return "Incorrect number of rows in matrix"

        return ''

    def _valid_format(self, table_json):
        """Format must be the expected version"""
        if table_json['format'] != self._format_version:
            return "Invalid format '%s', must be '%s'" % (table_json['format'],
                                                          self._format_version)
        else:
            return ''

    def _valid_type(self, table_json):
        """Table must be a known table type"""
        if table_json['type'].lower() not in self.TableTypes:
            return "Unknown BIOM type: %s" % table_json['type']
        else:
            return ''

    def _valid_generated_by(self, table_json):
        """Validate the generated_by field"""
        if not table_json['generated_by']:
            return "'generated_by' is not populated"
        if not isinstance(table_json['generated_by'], unicode):
            return "'generated_by' is not a string"

        return ''

    def _valid_nullable_id(self, table_json):
        """Validate the table id"""
        # this is nullable and don't actually care what is in here
        return ''

    def _valid_id(self, record):
        """Validate id for a row or column"""
        if not record['id']:
            return "'id' in %s appears empty" % record
        else:
            return ''

    def _valid_metadata(self, record):
        """Validate the metadata field for a row or column"""
        # this is nullable and don't actually care what is in here
        if record['metadata'] is None:
            return ''
        if isinstance(record['metadata'], dict):
            return ''

        return "metadata is neither null or an object"

    def _valid_rows(self, table_json):
        """Validate the 'rows' under 'table'"""
        required_keys = [('id', self._valid_id),
                         ('metadata', self._valid_metadata)]
        required_by_type = {}
        required_keys.extend(required_by_type.get(table_json['type'].lower(), []))

        for idx,row in enumerate(table_json['rows']):
            for key, method in required_keys:
                if key not in row:
                    return "ROW IDX %d MISSING '%s' FIELD" % (idx,key)

                result = method(row)
                if len(result) > 0:
                    return result
        return ''

    def _valid_columns(self, table_json):
        """Validate the 'columns' under 'table'"""
        required_keys = [('id', self._valid_id),
                         ('metadata', self._valid_metadata)]
        required_by_type = {} 
        required_keys.extend(required_by_type.get(table_json['type'].lower(), []))

        for idx, col in enumerate(table_json['columns']):
            for key, method in required_keys:
                if key not in col:
                    return "COL IDX %d MISSING '%s' FIELD" % (idx,key)

                result = method(col)
                if len(result) > 0:
                    return result
        return ''

    def _valid_data(self, table_json):
        """Validate the 'matrix' under 'table'"""
        if table_json['matrix_type'].lower() == 'sparse':
            return self._valid_sparse_data(table_json)
        elif table_json['matrix_type'].lower() == 'dense':
            return self._valid_dense_data(table_json)
        else:
            return "Unknown matrix type"

CommandConstructor = TableValidator
