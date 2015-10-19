# ----------------------------------------------------------------------------
# Copyright (c) 2011-2015, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .installation_informer import show_install_info
from .table_summarizer import summarize_table
from .metadata_adder import add_metadata
import biom.parse
import biom.util

__all__ = ['summarize_table', 'add_metadata', 'summarize_table']

def write_biom_table(table, format, filepath):
    """Write table in specified format to filepath"""

    if format not in ['hdf5', 'json', 'tsv']:
        raise ValueError("Unknown file format")

    if format == 'hdf5' and not biom.util.HAVE_H5PY:
        format = 'json'

    if format == 'json':
        with open(filepath, 'w') as f:
            f.write(table.to_json(biom.parse.generatedby()))
    elif format == 'tsv':
        with open(filepath, 'w') as f:
            f.write(table)
            f.write('\n')
    else:
        import h5py

        with h5py.File(filepath, 'w') as f:
            table.to_hdf5(f, biom.parse.generatedby())
