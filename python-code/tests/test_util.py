#!/usr/bin/env python

from os.path import abspath, dirname, exists

from biom.unit_test import TestCase, main
from biom.util import (natsort, _natsort_key, flatten, unzip,
                       get_biom_project_dir, parse_biom_config_files)

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Sandra Smit",
               "Zongzhi Liu", "Micah Hamady", "Daniel McDonald",
               "Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class UtilTests(TestCase):
    def test_natsort(self):
        """natsort should perform numeric comparisons on strings

        test pulled from QIIME (http://qiime.org)
        """
        # string with alpha and numerics sort correctly
        s = 'sample1 sample2 sample11 sample12'.split()
        self.assertEqual(natsort(s), 
          'sample1 sample2 sample11 sample12'.split())
        s.reverse()
        self.assertEqual(natsort(s), 
          'sample1 sample2 sample11 sample12'.split())
        self.assertEqual(natsort(list('cba321')),list('123abc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdba')),list('abcd'))

        # string of ints sort correctly
        self.assertEqual(natsort(['11','2','1','0']),
                               ['0','1','2','11'])

        # strings of floats sort correctly
        self.assertEqual(natsort(['1.11','1.12','1.00','0.009']),
                               ['0.009','1.00','1.11','1.12'])

        # string of ints sort correctly
        self.assertEqual(natsort([('11','A'),('2','B'),('1','C'),('0','D')]),
                            [('0','D'),('1','C'),('2','B'),('11','A')])

    def test_unzip(self):
        """unzip(items) should be the inverse of zip(*items)
        
        method pulled from PyCogent (http://pycogent.sourceforge.net)
        """
        chars = [list('abcde'), list('ghijk')]
        numbers = [[1,2,3,4,5], [0,0,0,0,0]]
        strings = [["abcde", "fghij", "klmno"], ['xxxxx'] * 3]
        empty = [[]]

        lists = [chars, numbers, strings]
        zipped = [zip(*i) for i in lists]
        unzipped = [unzip(i) for i in zipped]

        for u, l in zip(unzipped, lists):
            self.assertEqual(u, l)

    def test_flatten_no_change(self):
        """flatten should not change non-nested sequences (except to list)

        test pulled from PyCogent (http://pycogent.sourceforge.net)
        """
        self.assertEqual(flatten('abcdef'), list('abcdef')) #test identities
        self.assertEqual(flatten([]), []) #test empty sequence
        self.assertEqual(flatten(''), []) #test empty string

    def test_flatten(self):
        """flatten should remove one level of nesting from nested sequences

        test pulled from PyCogent (http://pycogent.sourceforge.net)
        """
        self.assertEqual(flatten(['aa', 'bb', 'cc']), list('aabbcc'))
        self.assertEqual(flatten([1,[2,3], [[4, [5]]]]), [1, 2, 3, [4,[5]]])

    def test_get_biom_project_dir(self):
        """Getting the biom project directory functions as expected.

        Test pulled from QIIME (http://qiime.org).
        """
        # Do an explicit check on whether the file system containing
        # the current file is case insensitive. This is in response
        # to SF bug #2945548, where this test would fail on certain
        # unusual circumstances on case-insensitive file systems
        # because the case of abspath(__file__) was inconsistent.
        # (If you don't believe this, set case_insensitive_filesystem
        # to False, and rename your top-level biom-format directory as
        # Biom-format on OS X. That should cause this test to fail as
        # actual will be path/to/Biom-format and expected will be
        # path/to/biom-format.) Note that we don't need to change anything
        # in the get_biom_project_dir() function as if the 
        # file system is case insenstive, the case of the returned
        # string is irrelevant.
        case_insensitive_filesystem = \
         exists(__file__.upper()) and exists(__file__.lower())

        actual = get_biom_project_dir()

        # I base the expected here off the imported location of
        # biom/util.py here, to handle cases where either the user
        # has biom-format in their PYTHONPATH, or when they've installed it
        # with setup.py.
        # If util.py moves this test will fail -- that
        # is what we want in this case, as the get_biom_project_dir()
        # function would need to be modified.
        import biom.util
        util_py_filepath = abspath(abspath(biom.util.__file__))
        expected = dirname(dirname(dirname(util_py_filepath)))

        if case_insensitive_filesystem:
            # Make both lowercase if the file system is case insensitive.
            actual = actual.lower()
            expected = expected.lower()
        self.assertEqual(actual,expected)

    def test_parse_biom_config_files(self):
        """parse_biom_config_files functions as expected.

        Test pulled from QIIME (http://qiime.org).
        """
        fake_file1 = ['key1\tval1', 'key2 val2']
        fake_file2 = ['key2\tval3']
        actual = parse_biom_config_files([fake_file1, fake_file2])
        expected = {'key1':'val1', 'key2':'val3'}
        self.assertEqual(actual, expected)

        # Looking up a nonexistent value returns None.
        self.assertEqual(actual['fake_key'], None)

        # Empty dict on empty input.
        self.assertEqual(parse_biom_config_files([]), {})


if __name__ == '__main__':
    main()
