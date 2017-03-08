#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import logging
import primer3
import starseqr_utils as su

logger = logging.getLogger("STAR-SEQR")


def runp3(seq_id, sequence, target=None):
    ''' can either use ":" in sequence, a base target or choose the middle of a read'''
    # get index of split to design targets
    if target:
        if (target + 30 > len(sequence)) or (target - 30 < 0):
            # target to close to border for primer generation
            return ()
        brk_target = [target - 10, 20]
    else:
        if ":" in sequence:
            mybrk = int(sequence.index(":"))
            sequence = sequence.replace(":", "")
            if (mybrk + 30 > len(sequence)) or (mybrk - 30 < 0):
                # target to close to border for primer generation
                return ()
            brk_target = [mybrk - 10, 20]
        else:
            if len(str(sequence)) < 75:
                return ()
            brk_target = [int(len(sequence) / 2), 1]
    # default values
    mydres = {
        'SEQUENCE_ID': seq_id,
        'SEQUENCE_TEMPLATE': sequence,
        # 'SEQUENCE_INCLUDED_REGION': [],
        'SEQUENCE_TARGET': brk_target  # start, len
    }
    mypres = {
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 27,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_OPT_GC_PERCENT': 50,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_LIBERAL_BASE': 0,
        'PRIMER_MIN_THREE_PRIME_DISTANCE': 0,  # can be 0=unique, or -1 multiple of same
        'PRIMER_LOWERCASE_MASKING': 0,  # 0 allows all case, 1, rejects if near 3'.
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PAIR_MAX_DIFF_TM': 6,
        'PRIMER_PRODUCT_SIZE_RANGE': [[75, 125], [125, 150],
                                      [150, 200], [200, 300]],
    }
    try:
        p3output = primer3.bindings.designPrimers(mydres, mypres)
        if (p3output['PRIMER_PAIR_NUM_RETURNED'] < 1):
            return ()
        else:
            p3res = parsep3(p3output)
    except (OSError):
        logger.error("Primer Design Failed- continuing", exc_info=True)
        return ()
    except (Exception):
        return ()
    return p3res


def parsep3(p3output):
    # logger.debug('Parsing Primer3 Results')
    Lprimer = str(p3output['PRIMER_LEFT_0_SEQUENCE'])
    # Ltm = round(p3output['PRIMER_LEFT_0_TM'])
    # Ltuple = p3output['PRIMER_LEFT_0']
    # Lstart, Llen = str(Ltuple[0]), str(Ltuple[1])
    Rprimer = str(p3output['PRIMER_RIGHT_0_SEQUENCE'])
    # Rtm = round(p3output['PRIMER_RIGHT_0_TM'])
    # Rtuple = p3output['PRIMER_RIGHT_0']
    # Rstart, Rlen = str(Rtuple[0]), str(Rtuple[1])
    # amplen = str(int(Rstart) - int(Lstart))
    return (Lprimer.upper(), Rprimer.upper())


def wrap_runp3(jxn, max_trx_fusion, chim_dir):
    '''
    use the breakpoint info stored in the fasta file header.
    Design primers using predicted fusion transcripts from expression
    '''
    clean_jxn = su.common.safe_jxn(jxn)
    fusionfa = os.path.join(chim_dir, 'transcripts-fusion-' + clean_jxn + '.fa')
    fusions_list = list(su.common.fasta_iter(fusionfa))  # list of tuples containing name, seq
    if len(fusions_list) > 0 and len(max_trx_fusion) > 0:
        for fusion_name, fus_seq in fusions_list:
            jxn, side, fusion_name, brk = fusion_name.split('|')
            brk = int(brk)
            if fusion_name in max_trx_fusion:
                return runp3(fusion_name, fus_seq.upper(), brk)
    else:
        return ()
