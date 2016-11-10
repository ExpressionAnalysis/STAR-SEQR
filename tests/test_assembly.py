import unittest
import sys
sys.path.insert(0, '../')
import starseqr_utils
import os
import shutil

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class AssemblyTestCase(unittest.TestCase):
    '''Tests assembly'''

    def setUp(self):
        jxn = 'chr19_580462_pos_chr12_120876181_pos_1_0'
        clean_jxn = str(jxn).replace(':', '_')
        clean_jxn = str(clean_jxn).replace('+', 'pos')
        clean_jxn = str(clean_jxn).replace('-', 'neg')
        jxn_dir = 'support' + '/' + clean_jxn + '/'
        assem_dir = jxn_dir + 'assem_jxn'
        shutil.rmtree(assem_dir, ignore_errors=True)

    def tearDown(self):
        jxn = 'chr19_580462_pos_chr12_120876181_pos_1_0'
        clean_jxn = str(jxn).replace(':', '_')
        clean_jxn = str(clean_jxn).replace('+', 'pos')
        clean_jxn = str(clean_jxn).replace('-', 'neg')
        jxn_dir = 'support' + '/' + clean_jxn + '/'
        assem_dir = jxn_dir + 'assem_jxn'
        shutil.rmtree(assem_dir, ignore_errors=True)

    def test_assembly_v1(self):
        """test velvet"""
        test_seq = 'TTCCTCCCCGAGCCCATGGGCACGGCCAACATCCAGCTCCACG:CTCGCATGTG'
        res1 = starseqr_utils.run_assembly.get_assembly_seq('chr19_580462_pos_chr12_120876181_pos_1_0', test_seq, 'velvet')
        expected_seq = 'CGATCCTCCGCTCCTGAGGCCCCCACAATGAAGCAGTCGGACGCGTCTCCCCAAGAAAGGGTGGACTCCGACGACCAGTGGGGAGAGTACTCCTGCGTCTTCCTCCCCGAGCCCATGGGCACG'
        assert(res1[0] == expected_seq)

    # def test_assembly_v2(self):
    #     """test spades--not working for this example"""
    #     test_seq = 'TTCCTCCCCGAGCCCATGGGCACGGCCAACATCCAGCTCCACG:CTCGCATGTG'
    #     res2 = starseqr_utils.run_assembly.get_assembly_seq('chr19_580462_pos_chr12_120876181_pos_1_0', test_seq, 'spades', )
    #     print(res2)
    #     expected_seq = ''
    #     assert(res2[0] == expected_seq)


if __name__ == '__main__':
    unittest.main()
