#!/usr/bin/env python

from sys import argv
from gzip import open as gzip_open
from biom.parse import parse_biom_table

if __name__ == '__main__':
    table = parse_biom_table(gzip_open(argv[1]))
    md = dict([(i, {'foo': 10}) for i in table.SampleIds])
    foo = table.addSampleMetadata(md)
