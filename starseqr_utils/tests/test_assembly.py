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
        jxn = 'chr1_45268528_pos_chr4_107152937_neg_0_2'
        clean_jxn = su.common.safe_jxn(jxn)
        jxn_dir = os.path.join('support', clean_jxn)
        assem_dir = os.path.join(jxn_dir, 'assem_pair')
        shutil.rmtree(assem_dir, ignore_errors=True)

    def tearDown(self):
        jxn = 'chr1_45268528_pos_chr4_107152937_neg_0_2'
        clean_jxn = su.common.safe_jxn(jxn)
        jxn_dir = os.path.join('support', clean_jxn)
        assem_dir = os.path.join(jxn_dir, 'assem_pair')
        shutil.rmtree(assem_dir, ignore_errors=True)
        os.remove(os.path.join(jxn_dir, 'assembly_log.txt'))

    def test_assembly_v1(self):
        """test velvet"""
        res1 = su.run_assembly.get_assembly_info('chr1:45268528:+:chr4:107152937:-:0:2', 'velvet')
        expected_res = ('TTATACAGGTATTTGGGAATAAAAGCAGACATACATGCATAAGCCAAGGCTTCATTATTGAAGTTTAGATATAGGAATGGAGCACAAAGTGAGTCAAGACATGACACAGCCCAGTGACCATACATCCGCCTCAGGGCCGTGGCCCTGTCTCAGCAGCACTTCTGGAGCCACATAGTTGGGGGTGCCACAGATGGTCTTCTTCCTCTGCTCCGGAGGCTCCAACCGGGCTG', '214', 'ENST00000372201.4_1--ENST00000510927.5_1,ENST00000372201.4_1--ENST00000394708.6_1,ENST00000372201.4_1--ENST00000273980.9_1,ENST00000372201.4_1--ENST00000467183.6_1,ENST00000372201.4_1--ENST00000394706.7_1,ENST00000372201.4_1--ENST00000432496.6_1,ENST00000372201.4_1--ENST00000361687.8_1,ENST00000465443.5_1--ENST00000510927.5_1,ENST00000465443.5_1--ENST00000394708.6_1,ENST00000465443.5_1--ENST00000273980.9_1,ENST00000465443.5_1--ENST00000467183.6_1,ENST00000465443.5_1--ENST00000394706.7_1,ENST00000465443.5_1--ENST00000432496.6_1,ENST00000465443.5_1--ENST00000361687.8_1,ENST00000476731.1_1--ENST00000510927.5_1,ENST00000476731.1_1--ENST00000394708.6_1,ENST00000476731.1_1--ENST00000273980.9_1,ENST00000476731.1_1--ENST00000467183.6_1,ENST00000476731.1_1--ENST00000394706.7_1,ENST00000476731.1_1--ENST00000432496.6_1,ENST00000476731.1_1--ENST00000361687.8_1')
        assert(res1 == expected_res)


if __name__ == '__main__':
    unittest.main()
