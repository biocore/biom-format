#!/bin/sh

# Author: Daniel McDonald
# Copyright 2012, BIOM Format
# This code is licensed under the GPL

SCRIPTS_DIR=$1
TABLES_DIR=$2
RESULTS_FP=$3

usage () {
    echo "usage: bench.sh <path_to_scripts> <path_to_tables> <path_to_results>"
    exit 0;
}

if [ -z $SCRIPTS_DIR ]
then
    usage ;
fi

if [ -z $TABLES_DIR ]
then
    usage ;
fi

if [ -e $RESULTS_FP ]
then
    rm $RESULTS_FP
fi

for f in `/bin/ls $SCRIPTS_DIR | grep "\.py$"`
do
    for t in `/bin/ls $TABLES_DIR`
    do
        /usr/bin/time -f "%M\t%e" python $SCRIPTS_DIR/$f $TABLES_DIR/$t 2> /tmp/foo
        if [ $? -ne 0 ]
        then
            echo "$f\t$t\tERROR" >> $RESULTS_FP
        else
            foo=`cat /tmp/foo`
            echo "$f\t$t\t$foo" >> $RESULTS_FP
        fi
    done
done
