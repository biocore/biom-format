#!/usr/bin/env python

"""Print information on the biom-format installation.

For more details, see http://biom-format.org
"""
from sys import platform, version as python_version, executable

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jose Clemente", "Greg Caporaso", 
               "Jai Ram Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", 
               "Sue Huse"]
__url__ = "http://biom-format.org"
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

try:
    from numpy import __version__ as numpy_lib_version
except ImportError:
    numpy_lib_version = "ERROR: Not installed - this is required! (This will also cause the BIOM library to not be importable.)"

try:
    from biom import __version__ as biom_lib_version
except ImportError:
    biom_lib_version = "ERROR: Can't find the BIOM library code (or numpy) - is it installed and in your $PYTHONPATH?"

try:
    from biom.exception import InvalidSparseBackendException
    from biom.table import SparseObj
except ImportError:
    SparseObj = "ERROR: Can't find the BIOM library code (or numpy) - is it installed and in your $PYTHONPATH?"
except InvalidSparseBackendException as e:
    SparseObj = "ERROR: %s" % e

def get_script_version():
    return __version__

def print_biom_config():
    system_info = [
     ("Platform", platform),
     ("Python/GCC version",python_version.replace('\n', ' ')),
     ("Python executable",executable)]
    
    max_len =  max([len(e[0]) for e in system_info])
    print "\nSystem information"
    print  "==================" 
    for v in system_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])

    version_info = [
     ("NumPy version", numpy_lib_version),
     ("biom-format library version", biom_lib_version),
     ("biom-format script version", get_script_version()),]
    
    max_len =  max([len(e[0]) for e in version_info])

    print "\nDependency versions"
    print  "===================" 
    for v in version_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])
    
    package_info = [
     ("SparseObj type", SparseObj)
    ]
    max_len =  max([len(e[0]) for e in package_info])
    
    print "\nbiom-format package information"
    print   "==============================="
    for v in package_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])
    print ""


if __name__ == '__main__':
    print_biom_config()
