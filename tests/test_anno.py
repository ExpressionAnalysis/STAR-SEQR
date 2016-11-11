import unittest
import os
import sys
import gzip
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import pandas as pd
sys.path.insert(0, '../')
import starseqr_utils


path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class AnnoTestCase(unittest.TestCase):
    """Tests annotation"""
    def find_resource(filename):
        packagedir = starseqr_utils.__path__[0]
        dirname = os.path.join(packagedir, 'resources')
        fullname = os.path.abspath(os.path.join(dirname, filename))
        return fullname

    ensgene = find_resource("ensGene.txt.gz")
    kg = gzip.open(ensgene)
    global gtree
    gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)

    def test_get_jxnside_anno_v1(self):
        """test get_jxnside_anno"""
        jxn_filt = pd.DataFrame({'name': 'chr2:42528535:+:chr2:29446451:-:0:0', 'jxn_reads': 'a1,a2,a3', 'jxn_counts': 3, 'spans': 5, 'spanreads': 's1,s2,s3,s4,s5', 'dist': 13082084, 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'}, index=pd.Series(0))
        jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'], jxn_filt['left_cdslen'], jxn_filt['left_exons'] = zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
        jxn_filt['right_symbol'], jxn_filt['right_annot'], jxn_filt['right_strand'], jxn_filt['right_cdslen'], jxn_filt['right_exons'] = zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_anno(x['name'], gtree, 2), axis=1))
        self.assertTrue(jxn_filt['left_symbol'].iloc[0]=='ENSG00000143924')
        self.assertTrue(jxn_filt['right_symbol'].iloc[0]=='ENSG00000171094')
        self.assertTrue(jxn_filt['left_strand'].iloc[0]=='+')
        self.assertTrue(jxn_filt['right_strand'].iloc[0]=='-')
        self.assertTrue(jxn_filt['left_cdslen'].iloc[0]==160596)
        self.assertTrue(jxn_filt['right_cdslen'].iloc[0]==727436)

    def test_get_jxnside_anno_v2(self):
        """test get_jxnside_anno.. Tests where no transcripts found"""
        jxn_filt = pd.DataFrame({'name': 'chr19:560462:+:chr8:560462:-:0:0', 'jxn_reads': 'a1,a2,a3', 'jxn_counts': 3, 'spans': 5, 'spanreads': 's1,s2,s3,s4,s5', 'dist': 13082084, 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'}, index=pd.Series(0))
        jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'], jxn_filt['left_cdslen'], jxn_filt['left_exons'] = zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
        self.assertTrue(jxn_filt['left_symbol'].iloc[0]=='NA')
        self.assertTrue(jxn_filt['left_strand'].iloc[0]=='NA')
        self.assertTrue(jxn_filt['left_cdslen'].iloc[0]=='NA')
        self.assertTrue(jxn_filt['left_exons'].iloc[0]==['NA'])

    def test_get_jxn_genes(self):
        """test_get_jxn_genes"""
        jxn_filt = pd.DataFrame({'name': 'chr2:42528535:+:chr2:29446451:-:0:0', 'jxn_reads': 'a1,a2,a3', 'jxn_counts': 3, 'spans': 5, 'spanreads': 's1,s2,s3,s4,s5', 'dist': 13082084, 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'}, index=pd.Series(0))
        jxn_filt['left_all'], jxn_filt['right_all'] = zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxn_genes(x['name'], gtree), axis=1))
        self.assertTrue(jxn_filt['left_all'].iloc[0]==['ENSG00000143924'])
        self.assertTrue(jxn_filt['right_all'].loc[0]==['ENSG00000171094'])

    def test_get_pos_genes(self):
        """test_get_pos_genes"""
        self.assertTrue(starseqr_utils.support_funcs_rna.get_pos_genes('chr2', 42528535, gtree)==['ENSG00000143924'])


if __name__ == '__main__':
    unittest.main()