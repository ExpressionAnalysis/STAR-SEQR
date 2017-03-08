import unittest
import os
import sys
import pandas as pd
sys.path.insert(0, '../')
import starseqr_utils as su

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class AnnoDBTestCase(unittest.TestCase):
    '''Tests annotation for db'''
    chimerdb3_path = su.common.find_resource('ChimerDB3.0_ChimerSeq.bedpe')
    fuca_path = su.common.find_resource('FusionCancer.bedpe')

    global chimerdb3
    chimerdb3 = pd.read_csv(chimerdb3_path, sep="\t", header=None, low_memory=False, engine='c')
    chimerdb3.columns = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Fusion_pair', 'Counts',
                         'strand1', 'strand2', 'Source', 'Cancer', 'Features']

    global fuca
    fuca = pd.read_csv(fuca_path, sep="\t", header=None, low_memory=False, engine='c')
    fuca.columns = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Database_ID', 'Rate',
                    'strand1', 'strand2', 'Methods', 'Cancers', 'Fusion_name']

    def test_annodb_v1(self):
        """CCDC6-RET fusion. This position has two overlapping breakpoints in chimerdb, 12 counts"""
        res1 = su.annotate_db.get_chimerdb('chr10:61665880:-:chr10:43612032:+:0:0', chimerdb3)
        print(res1)
        # (12, 'TCGA-PRADA,TCGA-TopHat-Fusion', 'LUAD,THCA', 'Oncogene,KB,Pub')
        assert(res1[0] == 12)
        assert("TCGA-PRADA" in res1[1].split(","))
        assert("THCA" in res1[2].split(","))

    def test_annodb_v1_2(self):
        """CCDC6-RET fusion. This position has two overlapping breakpoints in chimerdb, 12 counts"""
        res1 = su.annotate_db.get_fusioncancerdb('chr10:61665880:-:chr10:43612032:+:0:0', fuca)
        assert(res1 == (('', '', '')))

    def test_annodb_v2(self):
        """TMPRSS2-ERG Fusion. This position has two overlapping breakpoints in chimerdb, 300 counts, 3 methods"""
        res1 = su.annotate_db.get_chimerdb('chr21:42880007:-:chr21:39817543:+:0:0', chimerdb3)
        # (329, 'TCGA-FusionScan,TCGA-PRADA,TCGA-TopHat-Fusion', 'PRAD', 'KB,Pub')
        assert(res1[0] == 329)
        assert("TCGA-PRADA" in res1[1].split(","))
        assert("PRAD" in res1[2].split(","))

    def test_annodb_v3(self):
        """BCR-ABL Fusion. Found in both DBs"""
        res1 = su.annotate_db.get_chimerdb('chr22:23632600:-:chr9:133729450:+:0:0', chimerdb3)
        print(res1)
        # (3, 'TCGA-FusionScan,ChiTaRs', 'nan,LAML', 'Kinase,Oncogene,KB,Pub')
        assert(res1[0] == 3)
        assert("TCGA-FusionScan" in res1[1].split(","))
        assert("LAML" in res1[2].split(","))

    def test_annodb_v3_2(self):
        """BCR-ABL Fusion. Found in both DBs"""
        res1 = su.annotate_db.get_fusioncancerdb('chr22:23632600:-:chr9:133729450:+:0:0', fuca)
        # (0.85, 'tophat,chimerascan,SOAPfuse,FusionMap', 'Chronic_myelogenous_leukemia')
        assert(res1[0] == 0.85)
        assert("tophat" in res1[1].split(","))
        assert("Chronic_myelogenous_leukemia" in res1[2].split(","))


if __name__ == '__main__':
    unittest.main()
