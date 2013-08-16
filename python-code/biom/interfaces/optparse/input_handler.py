#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from biom.util import biom_open
from biom.parse import parse_biom_table

def load_biom_table(biom_fp):
    with biom_open(biom_fp, 'U') as table_f:
        return (parse_biom_table(table_f), biom_open(biom_fp, 'U'))

        