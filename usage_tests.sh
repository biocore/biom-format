#!/bin/bash

biom show-install-info
for table in examples/*hdf5.biom; do echo ${table}; biom validate-table -i ${table}; done
for table in examples/*table.biom; do echo ${table}; biom validate-table -i ${table}; done;
python biom/assets/exercise_api.py examples/rich_sparse_otu_table_hdf5.biom sample
python biom/assets/exercise_api.py examples/rich_sparse_otu_table_hdf5.biom observation
sh biom/assets/exercise_cli.sh
