#!/usr/bin/env python

from biom.util import biom_open
from biom.parse import parse_biom_table

def load_biom_table(biom_fp):
    with biom_open(biom_fp, 'U') as table_f:
        return parse_biom_table(table_f)