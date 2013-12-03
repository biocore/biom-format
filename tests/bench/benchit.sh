#!/bin/sh

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

PSTAT="import pstats; pstats.Stats('merged_t1_t2.stats').sort_stats('cumul').print_stats(15)"

# 0% overlap
OUT=bench_0per.txt

python make_and_profile_biom_table.py 100 200 10 100 200 10 0 100 5 0 100 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT
python make_and_profile_biom_table.py 1000 2000 10 1000 2000 10 0 1000 5 0 1000 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT
python make_and_profile_biom_table.py 10000 20000 10 10000 20000 10 0 10000 5 0 10000 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT

# 50% overlap
OUT=bench_50per.txt

python make_and_profile_biom_table.py 0 100 10 0 100 10 0 100 5 0 100 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT
python make_and_profile_biom_table.py 0 1000 10 0 1000 10 0 1000 5 0 1000 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT
python make_and_profile_biom_table.py 0 10000 10 0 10000 10 0 10000 5 0 10000 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT

# 100% overlap
OUT=bench_100per.txt

python make_and_profile_biom_table.py 0 100 5 0 100 5 0 100 5 0 100 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT
python make_and_profile_biom_table.py 0 1000 5 0 1000 5 0 1000 5 0 1000 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT
python make_and_profile_biom_table.py 0 10000 5 0 10000 5 0 10000 5 0 10000 5 >> $OUT
python -c "$PSTAT" | grep "function calls" >> $OUT
echo "" >> $OUT

