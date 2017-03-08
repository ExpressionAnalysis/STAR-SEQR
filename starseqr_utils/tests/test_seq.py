import os
import sys
sys.path.insert(0, '../')
import starseqr_utils as su
import unittest
import shutil
import pandas as pd


path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class ExonSeqTestCase(unittest.TestCase):
    '''Tests exon seq retrieval'''

    def setUp(self):
        jxn = 'chr1:10:+:chr2:20:+:0:0'
        clean_jxn = su.common.safe_jxn(jxn)
        jxn_dir = 'support' + '/' + clean_jxn + '/'
        shutil.rmtree(jxn_dir, ignore_errors=True)
        os.mkdir(jxn_dir)
        shutil.rmtree("res", ignore_errors=True)

    def tearDown(self):
        jxn = 'chr1:10:+:chr2:20:+:0:0'
        clean_jxn = su.common.safe_jxn(jxn)
        jxn_dir = 'support' + '/' + clean_jxn + '/'
        shutil.rmtree(jxn_dir, ignore_errors=True)
        shutil.rmtree("res", ignore_errors=True)

    def test_exons2seq_v1(self):
        """test_exons2seq + strand"""
        fa_path = 'test_data/ex1.fa'
        exons = [[[[1, 'chr1', 1, 10, '+', 'ENST1'], [2, 'chr1', 11, 20, '+', 'ENST1'], [3, 'chr1', 21, 30, '+', 'ENST1']],
                  [[1, 'chr1', 1, 10, '-', 'ENST2'], [2, 'chr1', 21, 30, '-', 'ENST2'], [3, 'chr1', 31, 40, '-', 'ENST2']]]]
        df = pd.DataFrame({'name': 'chr1:10:+:chr2:20:+:0:0', 'left_trx_exons': exons})
        df.apply(lambda x: su.core.exons2seq(fa_path, x['left_trx_exons'], x['name'], "leftv1", out_dir="res"), axis=1)

    def test_exons2seq2_v2(self):
        """test_exons2seq for fusions"""
        fa_path = 'test_data/ex1.fa'
        exons = [[[[1, 'chr1', 1, 10, '+', 'ENST1'], [2, 'chr1', 11, 20, '+', 'ENST1'], [3, 'chr1', 21, 30, '+', 'ENST1']],
                  [[1, 'chr1', 1, 10, '-', 'ENST2'], [2, 'chr1', 21, 30, '-', 'ENST2'], [3, 'chr1', 31, 40, '-', 'ENST2']]]]
        exons2 = [[[[1, 'chr2', 1, 10, '+', 'ENST3'], [2, 'chr2', 11, 20, '+', 'ENST3'], [3, 'chr2', 21, 30, '+', 'ENST3']],
                   [[1, 'chr2', 1, 10, '-', 'ENST4'], [2, 'chr2', 21, 30, '-', 'ENST4'], [3, 'chr2', 31, 40, '-', 'ENST4']]]]
        df = pd.DataFrame({'name': 'chr1:10:+:chr2:20:+:0:0', 'left_exons': exons, 'right_exons': exons2})
        df.apply(lambda x: su.core.exons2seq(fa_path, x['left_exons'], x['name'], "fusionv2", x['right_exons'], out_dir="res"), axis=1)

    def test_exons2seq_v3(self):
        """test_exons2seq: tests empty"""
        fa_path = 'test_data/ex1.fa'
        exons = ['NA']
        df = pd.DataFrame({'name': 'chr1:10:+:chr2:20:+:0:0', 'left_trx_exons': exons})
        df.apply(lambda x: su.core.exons2seq(fa_path, x['left_trx_exons'], x['name'], "leftv3", out_dir="res"), axis=1)


if __name__ == '__main__':
    unittest.main()
