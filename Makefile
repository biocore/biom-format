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
	flake8 biom setup.py
