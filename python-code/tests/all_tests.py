#!/usr/bin/env python
#
# suite of biom format unit tests.
# run suite by executing this file
# addapted from PyCogent alltests.py
#
import doctest, biom.unit_test as unittest, sys, os

__author__ = "Jose Carlos Clemente Litran"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Jose Carlos Clemente Litran", "Daniel McDonald", "Greg Caporaso"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.0.0-dev"
__maintainer__ = "Jose Carlos Clemente Litran"
__email__ = "jose.clemente@gmail.com"

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.
    
    __import__ only imports the top-level module.
    
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

def suite():
    modules_to_test = [
        'test_parse',
        'test_table',
        'test_sparsedict',
        'test_util',
        'test_unit_test',
        ]
    try:
        import biom._sparsemat
        modules_to_test.append('test_sparsemat')
    except ImportError:
        pass

    alltests = unittest.TestSuite()
    
    for module in modules_to_test:
        test = unittest.findTestCases(my_import(module))
        alltests.addTest(test)
    return alltests

class BoobyTrappedStream(object):
    def __init__(self, output):
        self.output = output

    def write(self, text):
        self.output.write(text)
        raise RuntimeError, "Output not allowed in tests"
        
    def flush(self):
        pass
        
    def isatty(self):
        return False

if __name__ == '__main__':
    if '--debug' in sys.argv:
        s = suite()
        s.debug()
    else:
        orig = sys.stdout
        if '--output-ok' in sys.argv:
            sys.argv.remove('--output-ok')
        else:
            sys.stdout = BoobyTrappedStream(orig)
        try:
            unittest.main(defaultTest='suite', argv=sys.argv)
        finally:
            sys.stdout = orig
