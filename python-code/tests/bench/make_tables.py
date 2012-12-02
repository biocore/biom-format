#!/usr/bin/env python

from math import ceil
from biom.table import table_factory, SparseObj
from biom.csmat import CSMat

def n_nonzero_items(n,m,p):
    """get the number of nonzero items"""
    return int(ceil(n * m * p))

N_DIRECTION = 1
def get_next_row_index(n, cur_row):
    """returns the next row"""
    global N_DIRECTION
    if (cur_row + N_DIRECTION) >= n:
        N_DIRECTION = -1
    elif (cur_row + N_DIRECTION) < 0:
        N_DIRECTION = 1

    return cur_row + N_DIRECTION

M_DIRECTION = 1
def get_next_col_index(m, cur_col):
    """returns the next col"""
    global M_DIRECTION
    if (cur_col + M_DIRECTION) >= m:
        M_DIRECTION = -1
    elif (cur_col + M_DIRECTION) < 0:
        M_DIRECTION = 1

    return cur_col + M_DIRECTION

if __name__ == '__main__':
    from sys import argv, exit
    
    if len(argv) != 5:
        print "usage: python make_tables.py n m p outdir"
        exit()

    n = int(argv[1])
    m = int(argv[2])
    p = float(argv[3])
    outdir = argv[4]

    row_ids = map(str, range(n))
    col_ids = map(str, range(m))

    obj = SparseObj(n,m)

    cur_row = -1; cur_col = -1
    rcv = {}
    
    if isinstance(obj, CSMat):
        for i in range(n_nonzero_items(n,m,p)):
            cur_row = get_next_row_index(n, cur_row)
            cur_col = get_next_col_index(m, cur_col)
            obj._coo_values.append(1)
            obj._coo_rows.append(cur_row)
            obj._coo_cols.append(cur_col)
    else:
        for i in range(n_nonzero_items(n,m,p)):
            cur_row = get_next_row_index(n, cur_row)
            cur_col = get_next_col_index(m, cur_col)
            rcv[(cur_row, cur_col)] = 1
        obj.update(rcv)

    table = table_factory(obj, col_ids, row_ids)

    f = open("%s/%dx%dx%0.3f_bench.biom" % (outdir,n,m,p), 'w')
    
    f.write(table.getBiomFormatJsonString('make_table-bench'))
    f.close()
        
