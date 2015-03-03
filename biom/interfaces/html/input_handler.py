#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import json
from biom.parse import MetadataMap, parse_biom_table

__author__ = "Evan Bolyen"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Evan Bolyen", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "BSD"
__maintainer__ = "Evan Bolyen"
__email__ = "ebolyen@gmail.com"


def load_biom_table(table_f):
    """Return a parsed BIOM table."""
    return parse_biom_table(table_f)


def load_biom_table_with_file_contents(biom_f):
    """Return a BIOM table and the original open filehandle as a tuple.

    Useful when additional computation needs to be performed on the file
    contents, such as an MD5 sum.

    WARNING: this function does not close the open filehandle that it returns.
    Users of this function are responsible for closing the filehandle when done
    using it!
    """
    table = parse_biom_table(biom_f)
    if hasattr(biom_f, 'seek'):
        biom_f.seek(0)
    return table, biom_f


def load_json_document(f):
    """Return a parsed JSON object."""
    return json.load(f)


def load_metadata(lines):
    """Parse a sample/observation metadata file, return a ``MetadataMap``.

    If ``lines`` is ``None``, this function will return ``None``.
    """
    if lines is not None:
        return MetadataMap.from_file(lines)

    return None
