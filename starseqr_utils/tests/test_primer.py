import unittest
import os
import sys
sys.path.insert(0, '../')
import starseqr_utils as su


path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)


class PrimerTestCase(unittest.TestCase):
    '''Tests primer design'''
    def test_runp3_v1(self):
        """test_runp3 with : in middle"""
        seq = "CATCAGCAGAATTCCCGCCATAGGTTAACTCTTTGACATTCTTATTACCACTGTGGGTCTCTTTGGGATCCCAGGGCACTGTAGGCAACTAGAACAATGTCCTTGGACTTGCAGAACTCCAGGAGTTTGCTCTGGTTGAGGTAAGGGTGACATTCCACCTGTAA:TTTTTCAAGACCAAAAGGTGTCACAATATATGAGCCAATAACCACATGCAAATTTAAACTGGAAAATAGAAAAGAAGATAAAGAAGGACATTTTCTATGGAATTCTCCTTTCTGGTGCCTAAGGAGCTCGGCCAAACAACCTTCAAGAGTTTCAGGAGGTGGTGATTAT"
        primer_res = su.run_primer3.runp3('test1', seq)
        assert(len(primer_res) == 2)
        assert(len(primer_res[0]) > 0)
        assert(len(primer_res[1]) > 0)

    def test_runp3_v1_1(self):
        """same as v1 in middle but this time use the location"""
        seq = "CATCAGCAGAATTCCCGCCATAGGTTAACTCTTTGACATTCTTATTACCACTGTGGGTCTCTTTGGGATCCCAGGGCACTGTAGGCAACTAGAACAATGTCCTTGGACTTGCAGAACTCCAGGAGTTTGCTCTGGTTGAGGTAAGGGTGACATTCCACCTGTAATTTTTCAAGACCAAAAGGTGTCACAATATATGAGCCAATAACCACATGCAAATTTAAACTGGAAAATAGAAAAGAAGATAAAGAAGGACATTTTCTATGGAATTCTCCTTTCTGGTGCCTAAGGAGCTCGGCCAAACAACCTTCAAGAGTTTCAGGAGGTGGTGATTAT"
        primer_res = su.run_primer3.runp3('test1', seq, 164)
        assert(len(primer_res) == 2)
        assert(len(primer_res[0]) > 0)
        assert(len(primer_res[1]) > 0)

    def test_runp3_v2(self):
        """test_runp3 with no target"""
        seq = ":CATCAGCAGAATTCCCGCCATAGGTTAACTCTTTGACATTCTTATTACCACTGTGGGTCTCTTTGGGATCCCAGGGCACTGTAGGCAACTAGAACAATGTCCTTGGACTTGCAGAACTCCAGGAGTTTGCTCTGGTTGAGGTAAGGGTGACATTCCACCTGTAATTTTTCAAGACCAAAAGGTGTCACAATATATGAGCCAATAACCACATGCAAATTTAAACTGGAAAATAGAAAAGAAGATAAAGAAGGACATTTTCTATGGAATTCTCCTTTCTGGTGCCTAAGGAGCTCGGCCAAACAACCTTCAAGAGTTTCAGGAGGTGGTGATTAT"
        primer_res = su.run_primer3.runp3('test1', seq)
        assert(primer_res == ())

    def test_runp3_v3(self):
        """test_runp3 with short seq"""
        seq = "CATCAGCAGAATTCCCGCCATAG"
        primer_res = su.run_primer3.runp3('test1', seq)
        assert(primer_res == ())

    def test_runp3_v4(self):
        """test_runp3 with small target"""
        seq = "CATCAGCA:GAATTCCCGCCATAGGTTAACTCTTTGACATTCTTATTACCACTGTGGGTCTCTTTGGGATCCCAGGGCACTGTAGGCAACTAGAACAATGTCCTTGGACTTGCAGAACTCCAGGAGTTTGCTCTGGTTGAGGTAAGGGTGACATTCCACCTGTAATTTTTCAAGACCAAAAGGTGTCACAATATATGAGCCAATAACCACATGCAAATTTAAACTGGAAAATAGAAAAGAAGATAAAGAAGGACATTTTCTATGGAATTCTCCTTTCTGGTGCCTAAGGAGCTCGGCCAAACAACCTTCAAGAGTTTCAGGAGGTGGTGATTAT"
        primer_res = su.run_primer3.runp3('test1', seq)
        assert(primer_res == ())

    def test_runp3_v5(self):
        """test_runp3 with no target"""
        seq = "CATCAGCAGAATTCCCGCCATAGGTTAACTCTTTGACATTCTTATTACCACTGTGGGTCTCTTTGGGATCCCAGGGCACTGTAGGCAACTAGAACAATGTCCTTGGACTTGCAGAACTCCAGGAGTTTGCTCTGGTTGAGGTAAGGGTGACATTCCACCTGTAATTTTTCAAGACCAAAAGGTGTCACAATATATGAGCCAATAACCACATGCAAATTTAAACTGGAAAATAGAAAAGAAGATAAAGAAGGACATTTTCTATGGAATTCTCCTTTCTGGTGCCTAAGGAGCTCGGCCAAACAACCTTCAAGAGTTTCAGGAGGTGGTGATTAT"
        primer_res = su.run_primer3.runp3('test1', seq)
        assert(len(primer_res) == 2)
        assert(len(primer_res[0]) > 0)
        assert(len(primer_res[1]) > 0)

    def test_wrap_runp3_v1(self):
        """test wrapper for runp3"""
        fusions = 'ENST00000372201.4_1--ENST00000510927.5_1|988'
        res1 = su.run_primer3.wrap_runp3('chr1:45268528:+:chr4:107152937:-:0:2', fusions, "test_data/chim_transcripts/")
        assert(len(res1) == 2)
        assert(len(res1[0]) > 0)
        assert(len(res1[1]) > 0)

    def test_wrap_runp3_v2(self):
        """test wrapper for runp3. Empty fusion list"""
        fusions = ''
        res1 = su.run_primer3.wrap_runp3('chr1:45268528:+:chr4:107152937:-:0:2', fusions, "test_data/chim_transcripts/")
        assert(res1 == ())


if __name__ == '__main__':
    unittest.main()
