import unittest
import os
import sys
#from nose.tools import set_trace; set_trace()
sys.path.insert(0, '../')
import starseqr_utils as su

path = os.path.dirname(__file__)
if path != '':
    os.chdir(path)

class TranslateTestCase(unittest.TestCase):
    '''Tests codon translation to peptide'''

    def test_translation(self):
        """test translation of sequence in default frame"""
        sequence = 'ATGCTGATGACTGAGTAAT'
	assert(su.translate_sequence.transcript_to_peptide(sequence,0) == "MLMTE*")

    def test_translation_rna(self):
        """test translation of sequence containing U in default frame"""
        sequence_rna = 'AUGCUGAUGACUGAGUAA'
	assert(su.translate_sequence.transcript_to_peptide(sequence_rna,0) == "MLMTE*")	

    def test_translation_frame(self):
        """test translation of sequence in specified frame"""
        sequence_frame = 'TAUGCUGAUGACUGAGUAA'
	assert(su.translate_sequence.transcript_to_peptide(sequence_frame,1) == "MLMTE*")	

    def test_translation_error(self):
        """test translation of when sequence contains an invalid codon"""
        sequence_error = 'TAUGCUNAUGACUGAGUAA'
	self.assertRaises(LookupError, su.translate_sequence.transcript_to_peptide(sequence_error,1))		

if __name__ == '__main__':
    unittest.main()
