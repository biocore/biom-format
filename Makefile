# ----------------------------------------------------------------------------
# Copyright (c) 2013--, biom-format development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(WITH_DOCTEST), TRUE)
	TEST_COMMAND = pytest --doctest-modules --doctest-glob='*.pyx'
else
	TEST_COMMAND = pytest
endif

.PHONY: doc lint test

test:
	$(TEST_COMMAND)
	sh usage_tests.sh

lint:
	flake8 biom setup.py --exclude=biom/tests/long_lines.py

doc:
	$(MAKE) -C doc clean html
