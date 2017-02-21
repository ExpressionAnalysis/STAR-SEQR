#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import logging
import numpy as np
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')


def find_unique_overhangs(reads_fq):
    '''get unique overhang reads per sam flag'''
    # +/+ = 401/321, +/- = 385/321, -/+ = 401/337, -/- = 385/337
    rfq_gen = su.common.FastqParser(reads_fq)
    res_seqs = set()
    for rfq in rfq_gen:
        # read_tag = int(float(re.split('_', rfq.header)[-1]))  # eventually break these down into separate strand
        res_seqs.add(rfq.sequence)

    # classify overhang lengths
    res_seq_len = [len(i) for i in res_seqs]
    all_min20 = np.sum(i > 20 for i in res_seq_len)  # how many of the overhangs are > 20
    all_min35 = np.sum(i > 35 for i in res_seq_len)  # how many of the overhangs are > 35
    return (len(res_seqs), all_min20, all_min35)


def get_diversity(jxn):
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = os.path.join('support', clean_jxn)
    overhang_fq = os.path.join(jxn_dir, 'overhang.fastq')
    overhang_res = find_unique_overhangs(overhang_fq)
    return overhang_res
