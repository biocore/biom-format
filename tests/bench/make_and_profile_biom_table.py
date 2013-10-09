#!/usr/bin/env python

from biom.table import SparseObj, SparseOTUTable
from numpy.random import randint, random

from sys import argv

t1lo,t1ho,t1so,t1ls,t1hs,t1ss,t2lo,t2ho,t2so,t2ls,t2hs,t2ss = map(int, argv[1:])

def make_sparse_fill(t,to,ts):
    for o in to:
        for s in ts:
            if random() < 0.05:
                t[o,s] = 1

def make_table(lowobs,highobs,ostep,lowsamp,highsamp,sstep):
    #table_obs = list(set(randint(lowobs,highobs,nobs)))
    #table_samp = list(set(randint(lowsamp,highsamp,nsamp)))
    table_obs = range(lowobs, highobs, ostep)
    table_samp = range(lowsamp, highsamp, sstep)
    table_data = SparseObj(len(table_obs), len(table_samp))
    make_sparse_fill(table_obs, table_samp, table_data)
    return SparseOTUTable(table_data, table_samp, table_obs)

t1 = make_table(t1lo,t1ho,t1so,t1ls,t1hs,t1ss)
t2 = make_table(t2lo,t2ho,t2so,t2ls,t2hs,t2ss)
print "Obs in t1: %d, Samp in t1: %d" % (len(t1.ObservationIds), len(t1.SampleIds))
print "Obs in t2: %d, Samp in t2: %d" % (len(t2.ObservationIds), len(t2.SampleIds))
print "Num obs overlap: %d, Num samp overlap: %d" % (len(set(t1.ObservationIds).intersection(set(t2.ObservationIds))), \
                                                     len(set(t1.SampleIds).intersection(set(t2.SampleIds))))

#f = open('t1.biom','w')
#f.write(t1.getBiomFormatJsonString('proftest'))
#f.close()
#f = open('t2.biom','w')
#f.write(t2.getBiomFormatJsonString('proftest'))
#f.close()

from profile import run

run("t1.merge(t2)", "merged_t1_t2.stats")

import pstats; pstats.Stats('merged_t1_t2.stats').strip_dirs().sort_stats('cumul').print_stats()
