#!/usr/bin/env python
"""Define BIOM exceptions"""

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2020, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2020, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"


class BiomException(Exception):
    pass


class BiomParseException(BiomException):
    pass


class TableException(BiomException):
    pass


class DisjointIDError(BiomException):
    pass


class UnknownAxisError(TableException):
    def __init__(self, axis):
        super().__init__()
        self.args = ("Unknown axis '%s'." % axis,)


class UnknownIDError(TableException):
    def __init__(self, missing_id, axis):
        super().__init__()
        self.args = ("The %s ID '%s' could not be found in the BIOM table." %
                     (axis, missing_id),)


class InvalidSparseBackendException(BiomException):
    pass
