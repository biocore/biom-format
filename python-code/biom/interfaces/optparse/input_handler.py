#!/usr/bin/env python

from biom.parse import parse_biom_table

def load_biom_table(biom_fp):
    return parse_biom_table(open(biom_fp,'U'))