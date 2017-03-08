import unittest
import os
import sys
import pandas as pd
sys.path.insert(0, '../../')
import starseqr_utils as su

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class HomologyGraphTestCase(unittest.TestCase):
    '''Tests homology between fusions'''

    def test_compare2jxns(self):
        jxn1 = "17:41179310:+:16:16341251:+:2:2"
        jxn2 = "17:41179310:+:16:18558623:-:2:2"
        res1 = su.homology_graph.get_pairwise_hom(jxn1, jxn2, "test_data/chim_transcripts/", "left")
        assert(res1 == 1.0)

    def test_homgroups(self):
        ''' Three homologous, 1 not homologous'''
        df_dict = {1: {'TPM_Fusion': 1.93,
                       'TPM_Left': 0.201,
                       'TPM_Right': 1.912,
                       'name': '17:41179310:+:16:16341251:+:2:2'},
                   2: {'TPM_Fusion': 1.53,
                       'TPM_Left': 0.201,
                       'TPM_Right': 13.303,
                       'name': '17:41179310:+:16:18558623:-:2:2'},
                   3: {'TPM_Fusion': 1.53,
                       'TPM_Left': 0.201,
                       'TPM_Right': 1.912,
                       'name': '17:41179310:+:16:14942755:+:2:2'},
                   4: {'TPM_Fusion': 1.53,
                       'TPM_Left': 0.201,
                       'TPM_Right': 13.303,
                       'name': '17:41179310:+:16:20558623:-:2:2'}}
        df = pd.DataFrame.from_dict(df_dict, 'index')
        res1 = su.homology_graph.prune_homology_graph(df, "test_data/chim_transcripts")
        assert(res1 == ['17:41179310:+:16:18558623:-:2:2', '17:41179310:+:16:20558623:-:2:2', '17:41179310:+:16:14942755:+:2:2'])


if __name__ == '__main__':
    unittest.main()
