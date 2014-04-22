#!/usr/bin/env python

from sys import argv
from gzip import open as gzip_open
from biom.parse import parse_biom_table

if __name__ == '__main__':
    table = parse_biom_table(gzip_open(argv[1]))

    md = dict([(s_id, {'FOO': i % 4})
              for i, s_id in enumerate(table.SampleIds)])
    table.addSampleMetadata(md)

    foo = table.binSamplesByMetadata(lambda x: x['FOO'])
