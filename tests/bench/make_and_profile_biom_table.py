#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from biom.table import SparseObj, SparseOTUTable
from numpy.random import random

from sys import argv

t1lo, t1ho, t1so, t1ls, t1hs, t1ss, t2lo, t2ho, t2so, t2ls, t2hs, t2ss = map(
    int, argv[1:])


def make_sparse_fill(t, to, ts):
    for o in to:
        for s in ts:
            if random() < 0.05:
                t[o, s] = 1


def make_table(lowobs, highobs, ostep, lowsamp, highsamp, sstep):
    table_obs = range(lowobs, highobs, ostep)
    table_samp = range(lowsamp, highsamp, sstep)
    table_data = SparseObj(len(table_obs), len(table_samp))
    make_sparse_fill(table_obs, table_samp, table_data)
    return SparseOTUTable(table_data, table_samp, table_obs)

t1 = make_table(t1lo, t1ho, t1so, t1ls, t1hs, t1ss)
t2 = make_table(t2lo, t2ho, t2so, t2ls, t2hs, t2ss)
print "Obs in t1: %d, Samp in t1: %d" % (len(t1.observation_ids),
                                         len(t1.sample_ids))
print "Obs in t2: %d, Samp in t2: %d" % (len(t2.observation_ids),
                                         len(t2.sample_ids))
print "Num obs overlap: %d, Num samp overlap: %d" % (
    len(set(t1.observation_ids).intersection(set(t2.observation_ids))),
    len(set(t1.sample_ids).intersection(set(t2.sample_ids)))
    )


from profile import run

run("t1.merge(t2)", "merged_t1_t2.stats")

import pstats
pstats.Stats(
    'merged_t1_t2.stats').strip_dirs(
).sort_stats(
    'cumul').print_stats(
)
