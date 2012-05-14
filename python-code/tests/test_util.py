#!/usr/bin/env python

from biom.unit_test import TestCase, main
from biom.util import natsort, _natsort_key, flatten, unzip

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

if __name__ == '__main__':
    main()