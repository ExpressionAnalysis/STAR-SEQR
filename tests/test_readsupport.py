import unittest
import os
import sys
sys.path.insert(0, '../')
import starseqr_utils

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class ReadSupportTestCase(unittest.TestCase):
    '''Tests read support'''

    def test_readsupport_v1(self):
        """PLSCR1 and PLSCR2 test case"""
        pass


if __name__ == '__main__':
    unittest.main()