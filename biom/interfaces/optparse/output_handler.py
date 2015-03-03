#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import exists
from pyqi.core.exception import IncompetentDeveloperError
from pyqi.core.interfaces.optparse.output_handler import write_list_of_strings
from biom.parse import generatedby
from biom.util import HAVE_H5PY

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "BSD"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


def write_subsetted_biom_table(result_key, data, option_value=None):
    """Write a string to a file"""
    if option_value is None:
        raise IncompetentDeveloperError("Cannot write output without a "
                                        "filepath.")

    if exists(option_value):
        raise IOError("Output path '%s' already exists." % option_value)

    table, fmt = data

    if fmt not in ['hdf5', 'json']:
        raise IncompetentDeveloperError("Unknown file format")

    if fmt == 'json':
        write_list_of_strings(result_key, table, option_value)
    else:
        if HAVE_H5PY:
            import h5py
        else:
            # This should never be raised here
            raise ImportError("h5py is not available, cannot write HDF5!")

        with h5py.File(option_value, 'w') as f:
            table.to_hdf5(f, generatedby())


def write_biom_table(result_key, data, option_value=None):
    """Write a string to a file"""
    if option_value is None:
        raise IncompetentDeveloperError("Cannot write output without a "
                                        "filepath.")

    if exists(option_value):
        raise IOError("Output path '%s' already exists." % option_value)

    table, fmt = data

    if fmt not in ['hdf5', 'json', 'tsv']:
        raise IncompetentDeveloperError("Unknown file format")

    if fmt == 'hdf5' and not HAVE_H5PY:
        fmt = 'json'

    if fmt == 'json':
        with open(option_value, 'w') as f:
            f.write(table.to_json(generatedby()))
    elif fmt == 'tsv':
        with open(option_value, 'w') as f:
            f.write(table)
            f.write('\n')
    else:
        import h5py

        with h5py.File(option_value, 'w') as f:
            table.to_hdf5(f, generatedby())
