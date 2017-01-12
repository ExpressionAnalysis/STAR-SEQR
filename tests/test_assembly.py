import unittest
import sys
sys.path.insert(0, '../')
import starseqr_utils as su
import os
import shutil

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class AssemblyTestCase(unittest.TestCase):
    '''Tests assembly'''

    def setUp(self):
        jxn = 'chr1_45268528_pos_chr4_107152937_neg_0_2]'
        clean_jxn = su.common.safe_jxn(jxn)
        jxn_dir = 'support' + '/' + clean_jxn + '/'
        assem_dir = jxn_dir + 'assem_jxn'
        shutil.rmtree(assem_dir, ignore_errors=True)

    def tearDown(self):
        jxn = 'chr1_45268528_pos_chr4_107152937_neg_0_2]'
        clean_jxn = su.common.safe_jxn(jxn)
        jxn_dir = 'support' + '/' + clean_jxn + '/'
        assem_dir = jxn_dir + 'assem_jxn'
        shutil.rmtree(assem_dir, ignore_errors=True)

    def test_assembly_v1(self):
        """test velvet"""
        # test_seq = 'TTCCTCCCCGAGCCCATGGGCACGGCCAACATCCAGCTCCACG:CTCGCATGTG'
        res1 = su.run_assembly.get_assembly_info('chr1:45268528:+:chr4:107152937:-:0:2', 'velvet')
        expected_res = ('CAGCCCGGTTGGAGCCTCCGGAGCAGAGGAAGAAGACCATCTGTGGCACCCCCAACTATGTGGCTCCAGAAGTGCTGCTGAGACAGGGCCACGGCCCTGAGGCGGATGTATGGTCACTGGGCTGTGTCATGTCTTGACTCACTTTGTGCTCCATTCCTATATCTAAACTTCAATAATGAAGCCTTGGCTTATGCATGTATGTCTGCTTTTATTCCCAAATACCTGTATAA', '214', 'ENST00000372201.4--ENST00000432496.6,ENST00000372201.4--ENST00000394708.6,ENST00000372201.4--ENST00000273980.9,ENST00000372201.4--ENST00000361687.8,ENST00000372201.4--ENST00000394706.7')
        assert(res1 == expected_res)


if __name__ == '__main__':
    unittest.main()
