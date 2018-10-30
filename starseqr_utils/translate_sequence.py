#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import logging
import sys

logger = logging.getLogger('STAR-SEQR')

def main():
    if len(sys.argv)>1:
     	sequence = sys.argv[1]
	frame = 0
        if len(sys.argv)>2:
	    try:
	        frame = int(sys.argv[2])
	    except Exception as err:
                print("Invalid frame. Defaulting to 0: " + str(err))
		print(transcript_to_peptide(sequence, frame))
            else:
                if frame > 2:
                    frame = 2
		    print("Fame cannot be higher than 2. Using 2")
		print(transcript_to_peptide(sequence, frame))
        else:
	    print(transcript_to_peptide(sequence, frame))

def _codon_to_AA(codon):
    _codon_table = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    codon = codon.upper().replace("U","T")
    try:
    	aa = _codon_table[codon]
    except KeyError as err:
	raise KeyError('Invalid codon included: %s' % codon) 
    else:
        if aa is None:
	    raise KeyError('Invalid codon included: %s' % codon) 
        else:
    	    return _codon_table[codon]


def transcript_to_peptide(sequence, frame=0):
    peptide = ""
    sequence = sequence[frame:]
    for codon in [sequence[i:i+3] for i in range(0, len(sequence), 3)]:
        if len(codon) == 3:
	    try:
	        amino_acid = _codon_to_AA(codon)
            except KeyError as err:
                logger.error("Exception:" + str(err))
                peptide = None
		break
            else:
                peptide = peptide + amino_acid
    return peptide

if __name__ == "__main__":
    main()
