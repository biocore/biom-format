#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import json
import sys
import h5py
from datetime import datetime
from operator import and_
from functools import reduce

import click
import numpy as np

from biom.cli import cli
from biom.util import biom_open, is_hdf5_file


@cli.command(name='validate-table')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input filepath to validate against the BIOM format'
                   ' specification')
@click.option('-f', '--format-version', default=None,
              help='The specific format version to validate against')
def validate_table(input_fp, format_version):
    """Validate a BIOM-formatted file.

    Test a file for adherence to the Biological Observation Matrix (BIOM)
    format specification. This specification is defined at
    http://biom-format.org

    Example usage:

    Validate the contents of table.biom for adherence to the BIOM format
    specification

    $ biom validate-table -i table.biom

    """
    valid, report = _validate_table(input_fp, format_version)
    click.echo("\n".join(report))
    if valid:
        # apparently silence is too quiet to be golden.
        click.echo("The input file is a valid BIOM-formatted file.")
        sys.exit(0)
    else:
        click.echo("The input file is not a valid BIOM-formatted file.")
        sys.exit(1)


def _validate_table(input_fp, format_version=None):
    result = TableValidator()(table=input_fp, format_version=format_version)
    return result['valid_table'], result['report_lines']


# Refactor in the future. Also need to address #664
class TableValidator:

    FormatURL = "http://biom-format.org"
    TableTypes = {
        'otu table',
        'pathway table',
        'function table',
        'ortholog table',
        'gene table',
        'metabolite table',
        'taxon table',
    }
    MatrixTypes = {'sparse', 'dense'}
    ElementTypes = {'int': int, 'str': str, 'float': float, 'unicode': str}
    HDF5FormatVersions = {(2, 0), (2, 0, 0), (2, 1), (2, 1, 0)}

    def run(self, **kwargs):
        is_json = not is_hdf5_file(kwargs['table'])

        if kwargs['format_version'] in [None, 'None']:
            if is_json:
                kwargs['format_version'] = '1.0.0'
            else:
                kwargs['format_version'] = '2.1'
        elif is_json:
            if kwargs['format_version'] != "1.0.0":
                raise ValueError("Only format 1.0.0 is valid for JSON")
        else:
            fmt_ver = [int(v) for v in kwargs['format_version'].split('.')]
            if tuple(fmt_ver) not in self.HDF5FormatVersions:
                raise ValueError("Unrecognized format version: %s" %
                                 kwargs['format_version'])

        with biom_open(kwargs['table']) as f:
            if is_json:
                try:
                    kwargs['table'] = json.load(f)
                except ValueError:
                    # if we hit this, we've already determined the file is
                    # not hdf5, so if it also is not json then it cannot be
                    # a valid biom-format file with the current formats.
                    raise ValueError("The provided table does not appear to "
                                     "be biom-format 1.0.0, 2.0.0 or 2.1.0.")
                return self._validate_json(**kwargs)
            else:
                kwargs['table'] = f

                if not isinstance(f, h5py.File):
                    print("Attempting to validate an HDF5 BIOM table, but the "
                          "table does not appear to be in HDF5 format!")
                    sys.exit(1)
                return self._validate_hdf5(**kwargs)

    def __call__(self, table, format_version=None):
        return self.run(table=table, format_version=format_version)

    def _validate_hdf5(self, **kwargs):
        table = kwargs['table']

        report_lines = []
        valid_table = True

        required_attrs = [
            ('format-url', self._valid_format_url),
            ('format-version', self._valid_hdf5_format_version),
            ('type', self._valid_type),
            ('shape', self._valid_shape),
            ('nnz', self._valid_nnz),
            ('generated-by', self._valid_generated_by),
            ('id', self._valid_nullable_id),
            ('creation-date', self._valid_creation_date)
        ]

        required_groups = ['observation', 'sample',
                           'observation/matrix', 'sample/matrix']

        required_datasets = ['observation/ids',
                             'observation/matrix/data',
                             'observation/matrix/indices',
                             'observation/matrix/indptr',
                             'sample/ids',
                             'sample/matrix/data',
                             'sample/matrix/indices',
                             'sample/matrix/indptr']

        for required_attr, attr_validator in required_attrs:
            if required_attr not in table.attrs:
                valid_table = False
                report_lines.append("Missing attribute: '%s'" % required_attr)
                continue

            status_msg = attr_validator(table)

            if len(status_msg) > 0:
                valid_table = False
                report_lines.append(status_msg)

        for group in required_groups:
            if group not in table:
                report_lines.append("Missing required '%s' group" % group)
                valid_table = False

        for dataset in required_datasets:
            if dataset not in table:
                report_lines.append("Missing required '%s' dataset" % dataset)
                valid_table = False

        if 'shape' in table.attrs:
            n_obs, n_samp = table.attrs['shape']
            obs_ids = table.get('observation/ids', None)
            samp_ids = table.get('sample/ids', None)

            if obs_ids is None:
                valid_table = False
                report_lines.append("observation/ids does not exist, cannot "
                                    "validate shape")

            if samp_ids is None:
                valid_table = False
                report_lines.append("sample/ids does not exist, cannot "
                                    "validate shape")

            if n_obs != len(obs_ids):
                valid_table = False
                report_lines.append("Number of observation IDs is not equal "
                                    "to the described shape")

            if n_samp != len(samp_ids):
                valid_table = False
                report_lines.append("Number of sample IDs is not equal "
                                    "to the described shape")
        else:
            report_lines.append("Missing 'shape' attribute")
            valid_table = False

        if 'format-version' in table.attrs:
            t_ver = '.'.join([str(v) for v in table.attrs['format-version']])
            if kwargs['format_version'] in ['2.0', '2.0.0']:
                if t_ver != '2.0':
                    error = "Table indicates it is version %s" % t_ver
                else:
                    report_lines.append("WARNING: 2.0 is not actively "
                                        "supported!")
                    error = self._valid_hdf5_metadata_v200(table)

                if error is not None:
                    report_lines.append(error)
            else:
                if t_ver != '2.1':
                    error = "Table indicates it is version %s" % t_ver
                else:
                    error = self._valid_hdf5_metadata_v210(table)

                if error is not None:
                    report_lines.append(error)

        return {'valid_table': valid_table, 'report_lines': report_lines}

    def _valid_hdf5_metadata_v200(self, table):
        no_md = np.array(["[]"])
        try:
            json.loads(table['observation'].get('metadata', no_md)[0])
        except ValueError:
            return ("Observation metadata do not appear to be formatted"
                    " correctly")

        try:
            json.loads(table['sample'].get('metadata', no_md)[0])
        except ValueError:
            return ("Sample metadata do not appear to be formatted"
                    " correctly")

    def _valid_hdf5_metadata_v210(self, table):
        if 'observation/metadata' not in table:
            return "Observation/metadata group is missing"
        if 'observation/group-metadata' not in table:
            return "Observation/group-metadata is missing"
        if 'sample/metadata' not in table:
            return "Sample/metadata group is missing"
        if 'sample/group-metadata' not in table:
            return "Sample/group-metadata is missing"

        n_obs_ids = len(table['observation/ids'])
        n_samp_ids = len(table['sample/ids'])
        for name, ds in table['observation/metadata'].items():
            if len(ds) != n_obs_ids:
                return "%s has %d entries, but expected %d" % (name, len(ds),
                                                               n_obs_ids)
        for name, ds in table['sample/metadata'].items():
            if len(ds) != n_samp_ids:
                return "%s has %d entries, but expected %d" % (name, len(ds),
                                                               n_samp_ids)

    def _validate_json(self, **kwargs):
        table_json = kwargs['table']

        # Need to make this an attribute so that we have this info during
        # validation.
        self._format_version = kwargs['format_version']

        report_lines = []
        valid_table = True

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

            status_msg = method(table_json)

            if len(status_msg) > 0:
                valid_table = False
                report_lines.append(status_msg)

        if 'shape' in table_json:
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

    def _json_or_hdf5_get(self, table, key):
        if hasattr(table, 'attrs'):
            item = table.attrs.get(key, None)
            if item is not None and isinstance(item, bytes):
                item = item.decode('utf8')
            return item
        else:
            return table.get(key, None)

    def _json_or_hdf5_key(self, table, key):
        if hasattr(table, 'attrs'):
            return key.replace('_', '-')
        else:
            return key

    def _is_int(self, x):
        """Return True if x is an int or numpy int"""
        return np.issubdtype(type(x), np.integer)

    def _valid_nnz(self, table):
        """Check if nnz seems correct"""
        if not self._is_int(table.attrs['nnz']):
            return "nnz is not an integer!"
        if table.attrs['nnz'] < 0:
            return "nnz is negative!"
        return ''

    def _valid_format_url(self, table):
        """Check if format_url is correct"""
        key = self._json_or_hdf5_key(table, 'format_url')
        value = self._json_or_hdf5_get(table, key)

        if value != self.FormatURL:
            return "Invalid '%s'" % key
        else:
            return ''

    def _valid_shape(self, table):
        """Matrix header is (int, int) representing the size of a 2D matrix"""
        a, b = self._json_or_hdf5_get(table, 'shape')

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

    def _valid_date(self, val):
        valid_times = ["%Y-%m-%d",
                       "%Y-%m-%dT%H:%M",
                       "%Y-%m-%dT%H:%M:%S",
                       "%Y-%m-%dT%H:%M:%S.%f"]
        if isinstance(val, bytes):
            val = val.decode('utf8')

        valid_time = False
        for fmt in valid_times:
            try:
                datetime.strptime(val, fmt)
                valid_time = True
                break
            except:  # noqa
                pass

        if valid_time:
            return ''
        else:
            return "Timestamp does not appear to be ISO 8601"

    def _valid_creation_date(self, table):
        """Verify datetime can be parsed

        Expects ISO 8601 datetime format (for example, 2011-12-19T19:00:00
                                          note that a 'T' separates the date
                                          and time)
        """
        return self._valid_date(table.attrs['creation-date'])

    def _valid_datetime(self, table):
        """Verify datetime can be parsed

        Expects ISO 8601 datetime format (for example, 2011-12-19T19:00:00
                                          note that a 'T' separates the date
                                          and time)
        """
        return self._valid_date(table['date'])

    def _valid_sparse_data(self, table_json):
        """All index positions must be integers and values are of dtype"""
        dtype = self.ElementTypes[table_json['matrix_element_type']]
        n_rows, n_cols = table_json['shape']
        n_rows -= 1  # adjust for 0-based index
        n_cols -= 1  # adjust for 0-based index

        for idx, coord in enumerate(table_json['data']):
            try:
                x, y, val = coord
            except:  # noqa
                return "Bad matrix entry idx %d: %s" % (idx, repr(coord))

            if not self._is_int(x) or not self._is_int(y):
                return "Bad x or y type at idx %d: %s" % (idx, repr(coord))

            if not isinstance(val, dtype):
                return "Bad value at idx %d: %s" % (idx, repr(coord))

            if x < 0 or x > n_rows:
                return "x out of bounds at idx %d: %s" % (idx, repr(coord))

            if y < 0 or y > n_cols:
                return "y out of bounds at idx %d: %s" % (idx, repr(coord))

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

    def _valid_hdf5_format_version(self, table):
        """Format must be the expected version"""
        ver = table.attrs['format-version']
        if tuple(ver) not in self.HDF5FormatVersions:
            return "Invalid format version '%s'" % str(ver)
        else:
            return ""

    def _valid_format(self, table_json):
        """Format must be the expected version"""
        formal = f"Biological Observation Matrix {self._format_version}"

        if table_json['format'] not in [formal, self._format_version]:
            return f"Invalid format '{table_json['format']}', must be '{self._format_version}'"  # noqa: E501
        else:
            return ''

    def _valid_type(self, table):
        """Table must be a known table type"""
        key = self._json_or_hdf5_key(table, 'type')
        value = self._json_or_hdf5_get(table, key)
        if value is None or value == "":
            return "Unknown table type, however that is likely okay."
        if value.lower() not in self.TableTypes:
            return "Unknown BIOM type: %s" % value
        else:
            return ''

    def _valid_generated_by(self, table):
        """Validate the generated_by field"""
        key = self._json_or_hdf5_key(table, 'generated_by')
        value = self._json_or_hdf5_get(table, key)
        if not value:
            return "'generated_by' is not populated"

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
        ttype = table_json.get('type')
        if ttype is None:
            ttype = ""
        required_keys.extend(
            required_by_type.get(ttype.lower(), []))

        for idx, row in enumerate(table_json['rows']):
            for key, method in required_keys:
                if key not in row:
                    return "ROW IDX %d MISSING '%s' FIELD" % (idx, key)

                result = method(row)
                if len(result) > 0:
                    return result
        return ''

    def _valid_columns(self, table_json):
        """Validate the 'columns' under 'table'"""
        required_keys = [('id', self._valid_id),
                         ('metadata', self._valid_metadata)]
        required_by_type = {}
        ttype = table_json.get('type')
        if ttype is None:
            ttype = ""
        required_keys.extend(
            required_by_type.get(ttype.lower(), []))

        for idx, col in enumerate(table_json['columns']):
            for key, method in required_keys:
                if key not in col:
                    return "COL IDX %d MISSING '%s' FIELD" % (idx, key)

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
