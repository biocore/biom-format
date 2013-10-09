#!/usr/bin/env python

from sys import argv
from gzip import open as gzip_open
from biom.parse import parse_biom_table
from random import shuffle

if __name__ == '__main__':
    table = parse_biom_table(gzip_open(argv[1]))

    ids = table.ObservationIds[:]
    shuffle(ids)
    to_keep = set(ids[:int(len(ids) / 2.0)])

    foo = table.filterObservations(lambda x,y,z: y in to_keep)
