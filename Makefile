# ----------------------------------------------------------------------------
# Copyright (c) 2013--, biom-format development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(WITH_DOCTEST), TRUE)
	TEST_COMMAND = python setup.py test -a --doctest-modules --doctest-glob='*.pyx'
else
	TEST_COMMAND = python setup.py test 
endif

test:
	$(TEST_COMMAND)
	biom show-install-info
	for table in examples/*hdf5.biom; do echo ${table}; biom validate-table -i ${table}; done
	for table in examples/*table.biom; do echo ${table}; biom validate-table -i ${table}; done;
	python biom/assets/exercise_api.py examples/rich_sparse_otu_table_hdf5.biom sample
	python biom/assets/exercise_api.py examples/rich_sparse_otu_table_hdf5.biom observation
	sh biom/assets/exercise_cli.sh

lint:
	flake8 biom setup.py

doc:
	$(MAKE) -C doc clean html
