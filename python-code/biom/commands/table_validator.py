#!/usr/bin/env python

from __future__ import division
import json
import dateutil.parser
from operator import and_
from pyqi.core.command import Command, Parameter, ParameterCollection
from pyqi.core.exception import CommandError

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Daniel McDonald", "Jose Clemente", "Greg Caporaso",
               "Jai Ram Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", "Sue Huse"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.2-dev"
__author__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class TableValidator(Command):
    BriefDescription = "Validate a BIOM-formatted file"
    LongDescription = ("Test a file for adherence to the Biological "
                       "Observation Matrix (BIOM) format specification. This "
                       "specification is defined at http://biom-format.org")

    Parameters = ParameterCollection([
        Parameter(Name='table_file', DataType=file,
                  Description='the input BIOM table file(-like) object',
                  Required=True),
        Parameter(Name='format_version', DataType=str,
                  Description='the specific format version to validate '
                  'against', Required=False,
                  Default='Biological Observation Matrix 1.0.0'),
        Parameter(Name='detailed_report', DataType=bool,
                  Description='include more details in the output report',
                  Required=False, Default=False)
    ])

    FormatURL = "http://biom-format.org"
    TableTypes = set(['otu table', 'pathway table', 'function table',
                      'ortholog table', 'gene table', 'metabolite table',
                      'taxon table'])
    MatrixTypes = set(['sparse', 'dense'])
    ElementTypes = {'int': int,'str': str,'float': float, 'unicode': unicode}

    def run(self, **kwargs):
        table_file = kwargs['table_file']
        # Need to make this an attribute so that we have this info during
        # validation.
        self._format_version = kwargs['format_version']
        detailed_report = kwargs['detailed_report']

        report_lines = []
        valid_table = True

        json_table = json.load(table_file)

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
            if key not in json_table:
                valid_table = False
                report_lines.append("Missing field: '%s'" % key)
                continue

            if detailed_report:
                report_lines.append("Validating '%s'..." % key)

            status_msg = method(json_table)

            if len(status_msg) > 0:
                valid_table = False
                report_lines.append(status_msg)

        if 'shape' in json_table:
            if detailed_report:
                report_lines.append("Validating 'shape' versus number of rows "
                                    "and columns...")

            if ('rows' in json_table and
                len(json_table['rows']) != json_table['shape'][0]):
                valid_table = False
                report_lines.append("Number of rows in 'rows' is not equal to "
                                    "'shape'")

            if ('columns' in json_table and
                len(json_table['columns']) != json_table['shape'][1]):
                valid_table = False
                report_lines.append("Number of columns in 'columns' is not "
                                    "equal to 'shape'")

        return {'valid_table': valid_table, 'report_lines': report_lines}

    def _is_int(self, x):
        """Return True if x is an int"""
        return isinstance(x, int)

    def _valid_format_url(self, json_table):
        """Check if format_url is correct"""
        if json_table['format_url'] != self.FormatURL:
            return "Invalid 'format_url'"
        else:
            return ''

    def _valid_shape(self, json_table):
        """A matrix header is (int, int) representing the size of a 2D matrix"""
        a,b = json_table['shape']

        if not (self._is_int(a) and self._is_int(b)):
            return "'shape' values do not appear to be integers"
        else:
            return ''

    def _valid_matrix_type(self, json_table):
        """Check if a valid matrix type exists"""
        if json_table['matrix_type'] not in self.MatrixTypes:
            return "Unknown 'matrix_type'"
        else:
            return ''

    def _valid_matrix_element_type(self, json_table):
        """Check if a valid element type exists"""
        if json_table['matrix_element_type'] not in self.ElementTypes:
            return "Unknown 'matrix_element_type'"
        else:
            return ''

    def _valid_datetime(self, json_table):
        """Verify datetime can be parsed

        Expects ISO 8601 datetime format (for example, 2011-12-19T19:00:00
                                          note that a 'T' separates the date 
                                          and time)
        """
        try:
            _ = dateutil.parser.parse(json_table['date'])
        except:
            return "Timestamp does not appear to be ISO 8601"
        else:
            return ''

    def _valid_sparse_data(self, json_table):
        """All index positions must be integers and values are of dtype"""
        dtype = self.ElementTypes[json_table['matrix_element_type']]
        n_rows, n_cols = json_table['shape']
        n_rows -= 1 # adjust for 0-based index
        n_cols -= 1 # adjust for 0-based index

        for idx, coord in enumerate(json_table['data']):
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

    def _valid_dense_data(self, json_table):
        """All elements must be of dtype and correspond to shape"""
        dtype = self.ElementTypes[json_table['matrix_element_type']]
        n_rows, n_cols = json_table['shape']

        for row in json_table['data']:
            if len(row) != n_cols:
                return "Incorrect number of cols: %s" % repr(row)

            if not reduce(and_, [isinstance(v, dtype) for v in row]):
                return "Bad datatype in row: %s" % repr(row)

        if len(json_table['data']) != n_rows:
            return "Incorrect number of rows in matrix"

        return ''

    def _valid_format(self, json_table):
        """Format must be the expected version"""
        if json_table['format'] != self._format_version:
            return "Invalid format '%s', must be '%s'" % (json_table['format'],
                                                          self._format_version)
        else:
            return ''

    def _valid_type(self, json_table):
        """Table must be a known table type"""
        if json_table['type'].lower() not in self.TableTypes:
            return "Unknown BIOM type: %s" % json_table['type']
        else:
            return ''

    def _valid_generated_by(self, json_table):
        """Validate the generated_by field"""
        if not json_table['generated_by']:
            return "'generated_by' is not populated"
        if not isinstance(json_table['generated_by'], unicode):
            return "'generated_by' is not a string"

        return ''

    def _valid_nullable_id(self, json_table):
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

    def _valid_rows(self, json_table):
        """Validate the 'rows' under 'table'"""
        required_keys = [('id', self._valid_id),
                         ('metadata', self._valid_metadata)]
        required_by_type = {}
        required_keys.extend(required_by_type.get(json_table['type'].lower(), []))

        for idx,row in enumerate(json_table['rows']):
            for key, method in required_keys:
                if key not in row:
                    return "ROW IDX %d MISSING '%s' FIELD" % (idx,key)

                result = method(row)
                if len(result) > 0:
                    return result
        return ''

    def _valid_columns(self, json_table):
        """Validate the 'columns' under 'table'"""
        required_keys = [('id', self._valid_id),
                         ('metadata', self._valid_metadata)]
        required_by_type = {} 
        required_keys.extend(required_by_type.get(json_table['type'].lower(), []))

        for idx, col in enumerate(json_table['columns']):
            for key, method in required_keys:
                if key not in col:
                    return "COL IDX %d MISSING '%s' FIELD" % (idx,key)

                result = method(col)
                if len(result) > 0:
                    return result
        return ''

    def _valid_data(self, json_table):
        """Validate the 'matrix' under 'table'"""
        if json_table['matrix_type'].lower() == 'sparse':
            return self._valid_sparse_data(json_table)
        elif json_table['matrix_type'].lower() == 'dense':
            return self._valid_dense_data(json_table)
        else:
            return "Unknown matrix type"

CommandConstructor = TableValidator
