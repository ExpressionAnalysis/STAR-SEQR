#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import logging
import re
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')


def find_unique_overhangs(reads_fq):
    '''get unique reads'''
    rfq_gen = su.common.FastqParser(reads_fq)
    # matches = {'321': 0, '337': 0, '385':0, '401': 0}
    res_left = set()
    res_right = set()
    for rfq in rfq_gen:
        read_tag = int(float(re.split('_', rfq.header)[-1]))
        if read_tag == 321 or read_tag == 337:
            res_left.add(rfq.sequence)
        else:
            res_right.add(rfq.sequence)
    return len(res_left), len(res_right)


def get_diversity(jxn):
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = 'support' + '/' + clean_jxn + '/'
    overhang_fq = jxn_dir + 'overhang.fastq'
    overhang_res = find_unique_overhangs(overhang_fq)
    return overhang_res
