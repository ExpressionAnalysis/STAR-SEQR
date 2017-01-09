#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import sys
import logging
import subprocess as sp
import re
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')


def do_velvet(assemdir, fastq, kmer, errlog, *args):
    ''' Run velvet with single or multiple fastqs '''
    logger.debug('*Running Velvet.')
    velveth_cmd = ['velveth', assemdir, str(kmer), '-short', '-fastq', fastq]
    for ar in args:
        velveth_cmd.extend(['-shortPaired', '-fastq', ar])
    vh_args = list(map(str, velveth_cmd))
    logger.debug('*velveth Command: ' + ' '.join(vh_args))
    try:
        p = sp.Popen(velveth_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(str(stdout))
        if stderr:
            errlog.write(str(stderr))
        if p.returncode != 0:
            logger.error('Error: velveth failed')
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('velveth Failed!', exc_info=True)
        sys.exit(1)
    velvetg_cmd = ['velvetg', assemdir, '-cov_cutoff', '2']
    vg_args = list(map(str, velvetg_cmd))
    logger.debug('*velvetg Command: ' + ' '.join(vg_args))
    try:
        p = sp.Popen(velvetg_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(str(stdout))
        if stderr:
            errlog.write(str(stderr))
        if p.returncode != 0:
            logger.error('Error: velvetg failed')
            sys.exit(1)
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('velvetg Failed!', exc_info=True)
        sys.exit(1)
    # extract header, sequence to generator
    records = su.common.fasta_iter(assemdir + '/contigs.fa')
    return records


# def do_spades(assemdir, pfastq, jxnfastq, errlog):
#     logger.debug('Running SPADES')
#     spades_cmd = ['spades.py', '--12', pfastq, '-s', jxnfastq,
#                   '-o', assemdir, '--phred-offset', 33, # --careful causing errors
#                   '-t', '1', '-m', '5',
#                   '--cov-cutoff', 'off']
#     spades_args = list(map(str, spades_cmd))
#     logger.debug('*SPADES Command: ' + ' '.join(spades_args))
#     errlog.write('*SPADES Command: ' + ' '.join(spades_args))
#     try:
#         p = sp.Popen(spades_args, stdout=sp.PIPE, stderr=sp.PIPE)
#         stdout, stderr = p.communicate()
#         if stdout:
#             errlog.write(stdout)
#         if stderr:
#             errlog.write(stderr)
#         if p.returncode != 0:
#             logger.error('Error: spades failed')
#     except (OSError) as o:
#         logger.error('Exception: ' + str(o))
#         logger.error('SPADES Failed!', exc_info=True)
#         sys.exit(1)
#     if (os.stat(os.path.realpath(assemdir + '/scaffolds.fasta')).st_size != 0):
#         records = fasta_iter(assemdir + '/scaffolds.fasta')
#         return records


def get_assembly_info(jxn, as_type):
    # clean jxn name to write to support folder made previous
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = 'support' + '/' + clean_jxn + '/'

    fusionfq = jxn_dir + 'transcripts_all_fusions.fa'
    fusions_list = list(su.common.fasta_iter(fusionfq))  # list of tuples containing name, seq

    pairfq = jxn_dir + 'paired.fastq'
    junctionfq = jxn_dir + 'junctions.fastq'

    # velvet
    assembly_list = []
    if as_type == 'velvet':
        errlog = open(jxn_dir + 'assembly_log.txt', 'w')
        assembly_list = list(do_velvet(jxn_dir + 'assem_pair', junctionfq, 17, errlog, pairfq))
        errlog.close()

    # elif as_type == 'spades':
    #     splog = open(jxn_dir + 'spades_log.txt', 'w')
    #     spades_all = do_spades(jxn_dir + 'spades', pairfq, junctionfq, splog)
    #     splog.close()

    # confirm assembly crosses breakpoint
    all_crossing = []
    all_seq = []
    all_len = []
    if len(assembly_list) > 0:
        for assembly in assembly_list:
            as_id, as_seq = assembly
            as_len = as_id.split('_')[3]
            # as_cov = int(float(as_id.split('_')[5]))
            all_seq.append(as_seq)
            all_len.append(as_len)
            # print(as_id, as_len, as_cov, as_seq)
            as_crossing_fusions = []
            if len(fusions_list) > 0:
                for fusion in fusions_list:
                    fusion_name, brk = fusion[0].split('|')
                    brk = int(brk)
                    fusion_seq = fusion[1][brk - 10:brk + 10].upper()
                    # print(fusion_name, brk, fusion_seq)
                    # regex solution from: http://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string
                    fusion_seq_re = re.compile('|'.join(fusion_seq[:i] + '.{0,2}' +
                                                        fusion_seq[i + 1:] for i in range(len(fusion_seq))))
                    if len(fusion_seq_re.findall(as_seq.upper())) or len(fusion_seq_re.findall(su.common.rc(as_seq).upper())) > 0:
                        as_crossing_fusions.append(fusion_name)
                all_crossing.extend(as_crossing_fusions)
    return all_seq, all_len, all_crossing
