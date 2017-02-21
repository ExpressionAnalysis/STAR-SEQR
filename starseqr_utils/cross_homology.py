#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import logging
import ssw
import numpy as np
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')


def get_ssw_scores(reads_fq, trxleft_fa, trxright_fa, mm_rate=2, mm_pct=.05):
    '''get sw score for left and right gene'''
    aligner = ssw.Aligner(gap_open=12, gap_extend=4)
    rfq_gen = su.common.FastqParser(reads_fq)
    matches = []
    for rfq in (x for _, x in zip(range(500), rfq_gen)):  # just do top 500 reads to maintain speed
        l_max = 0
        r_max = 0
        mod_len = min(rfq.seq_len - mm_rate, rfq.seq_len - rfq.seq_len * mm_pct)  # control error rates across a diversity of read len
        trxl_gen = su.common.fasta_iter(trxleft_fa)
        trxr_gen = su.common.fasta_iter(trxright_fa)
        for trxl_id, trxl_seq in trxl_gen:
            rfql_align = aligner.align(reference=rfq.sequence, query=trxl_seq)
            l_max = max(rfql_align.score, l_max)
        for trxr_id, trxr_seq in trxr_gen:
            rfqr_align = aligner.align(reference=rfq.sequence, query=trxr_seq)
            r_max = max(rfqr_align.score, r_max)
        read_norm = min(l_max, r_max) / mod_len  # pct identity matching
        matches.append(read_norm / 2)
    if len(matches) > 0:
        return(float("{0:.3f}".format(np.sum(i > .8 for i in matches) / len(matches))))
    else:
        return(0)


def get_cross_homology(jxn):
    paired_res = None
    overhang_res = None
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = os.path.join('support', clean_jxn)
    trans_left_fa = os.path.join(jxn_dir, 'transcripts_left.fa')
    trans_right_fa = os.path.join(jxn_dir, 'transcripts_right.fa')
    paired_fq = os.path.join(jxn_dir, 'span.fastq')
    overhang_fq = os.path.join(jxn_dir, 'overhang.fastq')
    paired_res = get_ssw_scores(paired_fq, trans_left_fa, trans_right_fa)
    overhang_res = get_ssw_scores(overhang_fq, trans_left_fa, trans_right_fa)
    return paired_res, overhang_res
