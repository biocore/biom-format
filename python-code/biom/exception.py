#!/usr/bin/env python

"""Define BIOM exceptions"""

class BiomException(Exception):
    pass

class BiomParseException(BiomException):
    pass

class TableException(BiomException):
    pass

class UnknownID(BiomException):
    pass
