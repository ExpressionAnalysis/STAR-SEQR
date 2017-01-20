import unittest
import os
import sys
import gzip
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import pandas as pd
sys.path.insert(0, '../')
import starseqr_utils as su


path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class AnnoTestCase(unittest.TestCase):
    """Tests annotation"""
    mygtf = os.path.realpath("test_data/test.gtf.gz")
    genepred_annot = os.path.splitext(mygtf)[0] + ".genePred"
    ucsc_annot = os.path.splitext(mygtf)[0] + ".UCSCTable.gz"
    su.gtf_convert.gtf_to_genepred(mygtf, genepred_annot)
    su.gtf_convert.genepred_to_UCSCtable(genepred_annot, ucsc_annot)
    su.gtf_convert.gtf_to_genepred(mygtf, genepred_annot)
    su.gtf_convert.genepred_to_UCSCtable(genepred_annot, ucsc_annot)
    kg = gzip.open(ucsc_annot)
    global gtree
    gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)

    def test_get_jxnside_anno_v1(self):
        """test get_jxnside_anno"""
        jxn_filt = pd.DataFrame({'name': 'chr1:871160:+:chr1:985950:-:0:0',
                                 'jxn_reads': 'a1,a2,a3',
                                 'jxn_counts': 3,
                                 'spans': 5,
                                 'spanreads': 's1,s2,s3,s4,s5',
                                 'dist': 114790,
                                 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'},
                                index=pd.Series(0))
        jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'], jxn_filt['left_cdslen'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
        jxn_filt['right_symbol'], jxn_filt['right_annot'], jxn_filt['right_strand'], jxn_filt['right_cdslen'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 2), axis=1))
        assert(jxn_filt['left_symbol'].iloc[0] == 'SAMD11')
        assert(jxn_filt['right_symbol'].iloc[0] == 'AGRN')
        assert(jxn_filt['left_strand'].iloc[0] == '+')
        assert(jxn_filt['right_strand'].iloc[0] == '+')
        assert(jxn_filt['left_cdslen'].iloc[0] == 9852)
        assert(jxn_filt['right_cdslen'].iloc[0] == 3395)

    def test_get_jxnside_anno_v2(self):
        """test get_jxnside_anno.. Tests where no transcripts found"""
        jxn_filt = pd.DataFrame({'name': 'chr19:560462:+:chr8:560462:-:0:0',
                                 'jxn_reads': 'a1,a2,a3',
                                 'jxn_counts': 3,
                                 'spans': 5,
                                 'spanreads': 's1,s2,s3,s4,s5',
                                 'dist': 13082084,
                                 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'},
                                index=pd.Series(0))
        jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'], jxn_filt['left_cdslen']= zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
        assert(jxn_filt['left_symbol'].iloc[0] == 'NA')
        assert(jxn_filt['left_strand'].iloc[0] == 'NA')
        assert(jxn_filt['left_cdslen'].iloc[0] == 'NA')

    def test_get_jxnside_anno_v3(self):
        """test get_jxnside_anno.. Tests non-coding transcripts for sorting"""
        jxn_filt = pd.DataFrame({'name': 'chr12:15704606:+:chr1:11918527:-:2:1',
                                 'jxn_reads': 'a1,a2,a3',
                                 'jxn_counts': 3,
                                 'spans': 5,
                                 'spanreads': 's1,s2,s3,s4,s5',
                                 'dist': 'NA',
                                 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'},
                                index=pd.Series(0))
        jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'], jxn_filt['left_cdslen'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
        assert(jxn_filt['left_symbol'].iloc[0] == 'NA')
        assert(jxn_filt['left_strand'].iloc[0] == 'NA')
        assert(jxn_filt['left_cdslen'].iloc[0] == 'NA')

    def test_get_jxn_genes(self):
        """test_get_jxn_genes"""
        jxn_filt = pd.DataFrame({'name': 'chr1:871160:+:chr1:985950:-:0:0',
                                 'jxn_reads': 'a1,a2,a3',
                                 'jxn_counts': 3,
                                 'spans': 5,
                                 'spanreads': 's1,s2,s3,s4,s5',
                                 'dist': 114790,
                                 'ann_format': 'Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length'},
                                index=pd.Series(0))
        jxn_filt['left_all'], jxn_filt['right_all'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxn_genes(x['name'], gtree), axis=1))
        assert(jxn_filt['left_all'].iloc[0] == ['SAMD11'])
        assert(jxn_filt['right_all'].loc[0] == ['AGRN'])

    def test_get_pos_genes(self):
        """test_get_pos_genes"""
        assert(su.annotate_sv.get_pos_genes('chr1', 871160, gtree) == ['SAMD11'])


if __name__ == '__main__':
    unittest.main()
