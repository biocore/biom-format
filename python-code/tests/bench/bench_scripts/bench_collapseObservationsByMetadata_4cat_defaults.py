#!/usr/bin/env python

from sys import argv
from gzip import open as gzip_open
from biom.parse import parse_biom_table

if __name__ == '__main__':
    table = parse_biom_table(gzip_open(argv[1]))
    
    md = [{'FOO': i % 4} for i in range(len(table.ObservationIds))]
    table.ObservationMetadata = md

    foo = table.collapseObservationsByMetadata(lambda x: x['FOO'])
