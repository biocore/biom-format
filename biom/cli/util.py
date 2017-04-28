# ----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

import biom.util
import biom.parse


def write_biom_table(table, fmt, filepath):
    """Write table in specified format to filepath"""

    if fmt not in ['hdf5', 'json', 'tsv']:
        raise ValueError("Unknown file format")

    if fmt == 'hdf5' and not biom.util.HAVE_H5PY:
        fmt = 'json'

    if fmt == 'json':
        with open(filepath, 'w') as f:
            f.write(table.to_json(biom.parse.generatedby()))
    elif fmt == 'tsv':
        with open(filepath, 'w') as f:
            f.write(table)
            f.write('\n')
    else:
        import h5py

        with h5py.File(filepath, 'w') as f:
            table.to_hdf5(f, biom.parse.generatedby())
