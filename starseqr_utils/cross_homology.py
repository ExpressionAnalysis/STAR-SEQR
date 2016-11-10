#!/usr/bin/env python

from __future__ import print_function
import logging
from itertools import groupby
import ssw
import gzip
import numpy as np

logger = logging.getLogger('STAR-SEQR')


def fasta_iter(fasta_name):
    ''' Given a fasta file, yield tuples of header, sequence'''
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
    for header in faiter:
        # drop the '>'
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = ''.join(s.strip() for s in faiter.next())
        # print(seq)
        yield header, seq
    fh.close()


class FastqRead(object):
    """Represents 1 read, or 4 lines of a casava 1.8-style fastq file."""
    def __init__(self, file_obj):
        header = file_obj.next().rstrip()
        assert header.startswith('@')
        self.header = header
        self.sequence = file_obj.next().rstrip()
        self.seq_len = len(self.sequence)
        file_obj.next()
        self.quality = file_obj.next().rstrip()
    def __str__(self):
        return "{0}\n{1}\n+\n{2}\n".format(self.header, self.sequence,
                                           self.quality, self.seq_len)


class FastqParser(object):
    """Parses a fastq file into FastqReads."""
    def __init__(self, filename, parse_headers=True):
        if filename.endswith('.gz'):
            self._file = gzip.open(filename, 'rb')
        else:
            self._file = open(filename, 'rU')
        self._line_index = 0
    def __iter__(self):
        return self
    def next(self):
        read = FastqRead(self._file)
        self._line_index += 4
        return read


def get_ssw_scores(reads_fq, trxleft_fa, trxright_fa):
    '''get sw score for left and right gene'''
    aligner = ssw.Aligner(gap_open=12, gap_extend=4)
    rfq_gen = FastqParser(reads_fq, parse_headers=False)
    matches = []
    for rfq in rfq_gen:
        l_max = 0
        r_max = 0
        trxl_gen = fasta_iter(trxleft_fa)
        trxr_gen = fasta_iter(trxright_fa)
        for trxl_id, trxl_seq in trxl_gen:
            rfql_align = aligner.align(reference=rfq.sequence, query=trxl_seq)
            l_max = max(rfql_align.score, l_max)
        for trxr_id, trxr_seq in trxr_gen:
            rfqr_align = aligner.align(reference=rfq.sequence, query=trxr_seq)
            r_max = max(rfqr_align.score, r_max)
        read_norm = min(l_max, r_max)
        matches.append(read_norm)
    if len(matches) > 0:
        return(int(np.percentile(matches, 75))) #UQ75 of ssw from cross-mapping
    else:
        return(None)


def get_cross_homology(jxn):
    # clean jxn name to write to support folder made previous
    paired_res = None
    overhang_res = None
    clean_jxn = str(jxn).replace(':', '_')
    clean_jxn = str(clean_jxn).replace('+', 'pos')
    clean_jxn = str(clean_jxn).replace('-', 'neg')
    jxn_dir = 'support' + '/' + clean_jxn + '/'

    trans_left_fa = jxn_dir + 'transcripts_left.fa'
    trans_right_fa = jxn_dir + 'transcripts_right.fa'
    paired_fq = jxn_dir + 'paired.fastq'
    overhang_fq = jxn_dir + 'overhang.fastq'
    paired_res = get_ssw_scores(paired_fq, trans_left_fa, trans_right_fa)
    overhang_res = get_ssw_scores(overhang_fq, trans_left_fa, trans_right_fa)
    return paired_res, overhang_res

