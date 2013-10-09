#!/usr/bin/env python

from sys import argv
from gzip import open as gzip_open
from biom.parse import parse_biom_table
from random import choice

if __name__ == '__main__':
    table = parse_biom_table(gzip_open(argv[1]))

    id_ = choice(table.SampleIds)
    foo = table.filterSamples(lambda x,y,z: y == id_)
