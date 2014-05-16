#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "BSD"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import exists
from pyqi.core.exception import IncompetentDeveloperError
from biom.parse import generatedby
from biom.util import HAVE_H5PY


def write_biom_table(result_key, data, option_value=None):
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
        with open(option_value, 'w') as f:
            f.write(table.to_json(generatedby()))
    else:
        if HAVE_H5PY:
            import h5py
        else:
            raise ImportError("h5py is not available, cannot write HDF5!")

        with h5py.File(option_value, 'w') as f:
            table.to_hdf5(f, generatedby())
