#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import logging
import subprocess as sp
from itertools import groupby
import re


logger = logging.getLogger('STAR-SEQR')


def fasta_iter(fasta_name):
    ''' Given a fasta file path, yield tuples of header, sequence '''
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


def do_velvet(assemdir, fastq, kmer, errlog, *args):
    ''' Run velvet with single or multiple fastqs '''
    logger.debug('*Running Velvet.')
    velveth_cmd = ['velveth', assemdir, str(kmer), '-short', '-fastq', fastq]
    for ar in args:
        velveth_cmd.extend(['-shortPaired', '-fastq', ar])
    vh_args = map(str, velveth_cmd)
    logger.debug('*velveth Command: ' + ' '.join(vh_args))
    try:
        p = sp.Popen(velveth_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error('Error: velveth failed')
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('velveth Failed!', exc_info=True)
        sys.exit(1)
    velvetg_cmd = ['velvetg', assemdir, '-cov_cutoff', '2']
    vg_args = map(str, velvetg_cmd)
    logger.debug('*velvetg Command: ' + ' '.join(vg_args))
    try:
        p = sp.Popen(velvetg_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error('Error: velvetg failed')
            sys.exit(1)
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('velvetg Failed!', exc_info=True)
        sys.exit(1)
    # extract header, sequence to generator
    records = fasta_iter(assemdir + '/contigs.fa')
    return records


def do_spades(assemdir, pfastq, jxnfastq, errlog):
    logger.debug('Running SPADES')
    spades_cmd = ['spades.py', '--12', pfastq, '-s', jxnfastq,
                  '-o', assemdir, '--phred-offset', 33, # --careful causing errors
                  '-t', '1', '-m', '5',
                  '--cov-cutoff', 'off']
    spades_args = map(str, spades_cmd)
    logger.debug('*SPADES Command: ' + ' '.join(spades_args))
    errlog.write('*SPADES Command: ' + ' '.join(spades_args))
    try:
        p = sp.Popen(spades_args, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error('Error: spades failed')
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('SPADES Failed!', exc_info=True)
        sys.exit(1)
    if (os.stat(os.path.realpath(assemdir + '/scaffolds.fasta')).st_size != 0):
        records = fasta_iter(assemdir + '/scaffolds.fasta')
        return records


def get_assembly_seq(jxn, jxn_seq, as_type):
    # TODO determine which transcript(s) have the best match to assembly by sw.
    # TODO Report transcript, crossing status, assembly quality
    # clean jxn name to write to support folder made previous
    clean_jxn = str(jxn).replace(':', '_')
    clean_jxn = str(clean_jxn).replace('+', 'pos')
    clean_jxn = str(clean_jxn).replace('-', 'neg')
    jxn_dir = 'support' + '/' + clean_jxn + '/'

    # get index of split to confine sequence
    predicted_seq = None
    if ":" in jxn_seq:
        mybrk = int(jxn_seq.index(":"))
        jxn_seq = jxn_seq.replace(":", "")
        predicted_seq = jxn_seq[mybrk-10:mybrk+10]

    pairfq = jxn_dir + 'paired.fastq'
    junctionfq = jxn_dir + 'junctions.fastq'
    final_seq = ''
    # velvet
    if as_type == 'velvet':
        errlog = open(jxn_dir + 'assembly_log.txt', 'w')
        velvet_all = do_velvet(jxn_dir + 'assem_pair', junctionfq, 17, errlog, pairfq)
        errlog.close()
        if velvet_all:
            for contig in velvet_all:
                final_id, final_seq = contig
                break
    elif as_type == 'spades':
        splog = open(jxn_dir + 'spades_log.txt', 'w')
        spades_all = do_spades(jxn_dir + 'spades', pairfq, junctionfq, splog)
        splog.close()
        if spades_all:
            for contig in spades_all:
                final_id, final_seq = contig
                break

    # confirm assembly crosses breakpoint
    found_status = 0
    if predicted_seq:
        predicted_seq = predicted_seq.upper()
        # regex solution from: http://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string
        predicted_seq_re=re.compile('|'.join(predicted_seq[:i]+'.{0,2}'+
                                             predicted_seq[i+1:] for i in range(len(predicted_seq))))
        if len(predicted_seq_re.findall(final_seq.upper())) > 0:
            found_status = 1
    return final_seq, found_status
