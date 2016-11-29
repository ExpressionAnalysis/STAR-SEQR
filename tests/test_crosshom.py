import unittest
import os
import sys
sys.path.insert(0, '../')
import starseqr_utils as su

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class CrossHomologyTestCase(unittest.TestCase):
    '''Tests cross homology'''

    def test_crosshom_v1(self):
        """test homology paired.fastq TUBA1B--TUBA1A test case"""
        res1 = su.cross_homology.get_cross_homology('chr12:49521770:-:chr12:49578857:-:6:0')
        assert(res1 == (91, 37))

    def test_crosshom_v2(self):
        """test homology paired.fastq no reads in fastq test case"""
        res1 = su.cross_homology.get_cross_homology('chr2:131389025:+:chr2:131996848:+:0:4')
        assert(res1 == (None, None))


if __name__ == '__main__':
    unittest.main()
