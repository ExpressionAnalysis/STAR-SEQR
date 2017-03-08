import unittest
import os
import sys
import gzip
sys.path.insert(0, '../')
import starseqr_utils as su

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class GTFTestCase(unittest.TestCase):
    '''Tests gtf conversion'''

    def test_gtf2genepred_v1(self):
        """test gtf to genepred"""
        mygtf = os.path.realpath("test_data/test.gtf.gz")
        genepred_annot = os.path.splitext(mygtf)[0] + ".genePred"
        os.remove(genepred_annot) if os.path.exists(genepred_annot) else None
        su.gtf_convert.gtf_to_genepred(mygtf, genepred_annot)
        assert os.path.exists(genepred_annot) == 1
        with open(genepred_annot) as myfile:
            count = sum(1 for line in myfile)
        assert(count) == 76

    def test_genepred2ucsc_v1(self):
        """test genepred to ucsc"""
        mygtf = os.path.realpath("test_data/test.gtf.gz")
        genepred_annot = os.path.splitext(mygtf)[0] + ".genePred"
        ucsc_annot = os.path.splitext(mygtf)[0] + ".UCSCTable.gz"
        os.remove(genepred_annot) if os.path.exists(genepred_annot) else None
        os.remove(ucsc_annot) if os.path.exists(ucsc_annot) else None
        su.gtf_convert.gtf_to_genepred(mygtf, genepred_annot)
        su.gtf_convert.genepred_to_UCSCtable(genepred_annot, ucsc_annot)
        assert os.path.exists(ucsc_annot) == 1
        with gzip.open(ucsc_annot, 'rb') as myfile:
            count = sum(1 for line in myfile)
        assert(count) == 76

    def test_genepred2ucsc_v2(self):
        """test genepred to ucsc. test for re-use"""
        mygtf = os.path.realpath("test_data/test.gtf.gz")
        genepred_annot = os.path.splitext(mygtf)[0] + ".genePred"
        ucsc_annot = os.path.splitext(mygtf)[0] + ".UCSCTable.gz"
        su.gtf_convert.gtf_to_genepred(mygtf, genepred_annot)
        su.gtf_convert.genepred_to_UCSCtable(genepred_annot, ucsc_annot)
        assert os.path.exists(ucsc_annot) == 1
        with gzip.open(ucsc_annot, 'rb') as myfile:
            count = sum(1 for line in myfile)
        assert(count) == 76


if __name__ == '__main__':
    unittest.main()
