#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from biom.util import biom_open
from biom.parse import parse_biom_table

def load_biom_table(biom_fp):
    """Return a parsed BIOM table."""
    with biom_open(biom_fp, 'U') as table_f:
        return parse_biom_table(table_f)

def load_biom_table_with_file_contents(biom_fp):
    """Return a BIOM table and the original open filehandle as a tuple.

    Useful when additional computation needs to be performed on the file
    contents, such as an MD5 sum.

    WARNING: this function does not close the open filehandle that it returns.
    Users of this function are responsible for closing the filehandle when done
    using it!
    """
    biom_f = biom_open(biom_fp, 'U')
    table = parse_biom_table(biom_f)
    biom_f.seek(0)
    return table, biom_f
