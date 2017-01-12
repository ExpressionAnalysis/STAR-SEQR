#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import six
import os
import sys
import re
import time
import logging
import subprocess as sp
import pysam  # requires 0.9.0
import collections
import numpy as np
import starseqr_utils as su

try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans

logger = logging.getLogger('STAR-SEQR')


def find_support_reads(jxn, bam, side, tx, sub_reads, gtree):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1 = int(pos1)
    pos2 = int(pos2)

    if side == 2:
        chrom1, chrom2 = chrom2, chrom1
        pos1, pos2 = pos2, pos1
        flipstr = maketrans('-+', '+-')
        str1, str2 = str2.translate(flipstr), str1.translate(flipstr)
        repleft, repright = repright, repleft

    maxhom = max(int(repright), int(repleft))
    if str1 == '+':
        pos1left = pos1 - 100000
        pos1right = pos1 + maxhom + 10
    if str2 == '-':
        pos2left = pos2 - 100000
        pos2right = pos2 + maxhom + 10
    if str1 == '-':
        pos1left = pos1 - maxhom - 10
        pos1right = pos1 + 100000
    if str2 == '+':
        pos2left = pos2 - maxhom - 10
        pos2right = pos2 + 100000

    # define all combinations upfront for recordkeeping
    retDict = {}
    lmode = ['spanleft', 'spanright', 'jxnleft', 'jxnright', 'hangleft', 'hangright']
    lstrand = ['for', 'rev']
    lorient = ['first', 'second']
    for x in lmode:
        retDict[x] = {}
        for y in lstrand:
            for z in lorient:
                retDict[x][(y, z)] = {}

    for read in bam:
        mymode = ''
        # determine if span, jxn read or chimeric overhang
        if (read.next_reference_name == chrom2 and read.reference_name == chrom1 and
                not (read.flag & 2) and not (read.flag & 256)):
            if ((int(read.reference_start) > pos1left and
                 int(read.reference_start) < pos1right) and
                (int(read.next_reference_start) > pos2left and
                 int(read.next_reference_start) < pos2right)):
                mymode = 'span'
        elif (read.next_reference_name == read.reference_name and
                (read.flag & 2) and not (read.flag & 256)):
            if ((int(read.reference_start) > pos1left and
                 int(read.reference_start) < pos1right) or
                (int(read.next_reference_start) > pos2left and
                 int(read.next_reference_start) < pos2right)):
                mymode = 'jxn'
        elif (read.next_reference_name == chrom2 and read.reference_name == chrom1 and
                not (read.flag & 2) and (read.flag & 256)):
            if ((int(read.reference_start) > pos1left and
                 int(read.reference_start) < pos1right) or
                (int(read.next_reference_start) > pos2left and
                 int(read.next_reference_start) < pos2right)):
                mymode = 'hang'

        if mymode:
            myside = 'left' if side == 1 else 'right'
            mytype = mymode + myside
            mystrand = 'for' if (read.flag & 16) else 'rev'
            myorient = 'first' if (read.flag & 64) else 'second'
            is_subset = 'SUB' if read.query_name in sub_reads else 'ALL'
            if is_subset == 'SUB':
                read_tx = su.annotate_sv.get_pos_genes(read.reference_name, read.reference_start, gtree)
                read_tx2 = su.annotate_sv.get_pos_genes(read.next_reference_name, read.next_reference_start, gtree)
                if set(tx).intersection(set(read_tx)) and set(tx).intersection(set(read_tx2)):
                    retDict[mytype][(mystrand, myorient)][read.query_name] = (is_subset,
                                                                              read.get_tag('AS'),
                                                                              read.get_tag('nM'),
                                                                              read.query_alignment_length,
                                                                              int(np.mean(read.query_alignment_qualities)))
    return retDict


def get_reads_from_bam(bam_file, jxn, tx, s_reads, gtree):
    '''
    Fetch reads from each direction of jxn and get an accounting of each if unique or duplicate.
    '''
    logger.debug('Extracting supporting read information for ' + jxn)
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)

    if not os.path.exists(bam_file):
        logger.error('BAM file could not be found: ' + bam_file)
        sys.exit()
    bamObject = pysam.Samfile(bam_file, 'rb')
    if not bamObject.check_index():
        logger.error('Index file could not be found for: ' + bam_file)
        sys.exit()

    # account for genome boundaries
    len1 = six.next((item for item in bamObject.header['SQ'] if item['SN'] == chrom1))['LN']
    len2 = six.next((item for item in bamObject.header['SQ'] if item['SN'] == chrom2))['LN']
    dist1_less = max(i for i in [int(pos1) - 100000, 0] if i >= 0)
    dist1_plus = min(i for i in [int(pos1) + 100000, len1] if i >= 0)
    dist2_less = max(i for i in [int(pos2) - 100000, 0] if i >= 0)
    dist2_plus = min(i for i in [int(pos2) + 100000, len2] if i >= 0)

    results = {}
    # support left
    logger.debug('Extracting left paired reads for ' + jxn)
    bam = bamObject.fetch(chrom1, dist1_less, dist1_plus)
    support_for = find_support_reads(jxn, bam, 1, tx, s_reads, gtree)
    results['spanleft'] = support_for['spanleft']
    results['jxnleft'] = support_for['jxnleft']
    results['hangleft'] = support_for['hangleft']
    # support right
    logger.debug('Extracting right paired reads for ' + jxn)
    bam = bamObject.fetch(chrom2, dist2_less, dist2_plus)
    support_rev = find_support_reads(jxn, bam, 2, tx, s_reads, gtree)
    results['spanright'] = support_rev['spanright']
    results['jxnright'] = support_rev['jxnright']
    results['hangright'] = support_rev['hangright']
    bamObject.close()
    logger.debug('Finished extracting supporting read information for ' + jxn)
    return results


def subset_bam_by_reads(bam, out_bam, read_ids, jxn):
    # bamfilternames is sensitive to duplicate read names in the list. Make sure they are unique.
    logger.debug('Subset bam with supporting reads to ' + out_bam)
    names = 'names=' + read_ids
    index = 'index=1'
    indexfilename = 'indexfilename=' + out_bam + '.bai'
    tmpfile = 'tmpfile=' + out_bam + '.tmp'
    stdinfile = bam
    bamf = open(out_bam, 'wb')
    try:
        with open(stdinfile, 'rb') as f:
            retcode = sp.call(['bamfilternames', names, index, indexfilename,
                               tmpfile], stdin=f, stdout=bamf, stderr=sp.PIPE)
        bamf.close()
        if retcode != 0:
            logger.error('bamfilternames failed on:' + jxn)
    except OSError as o:
        logger.error('bamfilternames Failed', exc_info=True)
        logger.error('Exception: ' + str(o))
        bamf.close()
        os._exit(1)


def bam2fastq(jxn_dir, in_bam, junctionfq, pairfq, overhangfq):
    '''
    STAR writes 3 lines in the SAM file for chimeric reads.
    Custom convert bam2fastq where pair1/2 go into paired.fastq
    and chimeric only reads go into junctions.fq.
    Chimeric reads need all seqs/quals while
    paired must be soft-clipped seqs/quals.
    '''
    logger.debug('Converting bam to fastq for ' + in_bam)
    if not os.path.exists(in_bam):
        logger.error('BAM file could not be found: ' + in_bam)
        sys.exit(1)
    pairfqfh = open(pairfq, 'w')
    junctionfqfh = open(junctionfq, 'w')
    overhangfqfh = open(overhangfq, 'w')
    try:
        bam_sort = in_bam[:-4] + '.nsorted'
        pysam.sort('-n', in_bam, '-o', bam_sort + '.bam')
        if not os.path.exists(bam_sort + '.bam'):
            logger.error('name sorted BAM file could not be found: ' + bam_sort + '.bam')
            sys.exit(1)
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('pysam Failed!', exc_info=True)
        sys.exit(1)
    bamObject = pysam.Samfile(bam_sort + '.bam', 'rb')
    for read in bamObject.fetch(until_eof=True):
        if not read.flag & 256:
            # print(read.query_alignment_sequence)
            if read.is_read1:
                orient = 1
                if read.is_reverse:
                    seq = su.common.rc(read.query_alignment_sequence)
                    quals = read.query_alignment_qualities[::-1]
                else:
                    seq = read.query_alignment_sequence
                    quals = read.query_alignment_qualities
            else:
                orient = 2
                if read.is_reverse:
                    seq = su.common.rc(read.query_alignment_sequence)
                    quals = read.query_alignment_qualities[::-1]
                else:
                    seq = read.query_alignment_sequence
                    quals = read.query_alignment_qualities
            pairfqfh.write('@' + read.query_name + '/' + str(orient) + '\n')
            pairfqfh.write(seq + '\n')
            pairfqfh.write('+' + '\n')
            pairfqfh.write(''.join(
                list(map(chr, [x + 33 for x in quals]))) + '\n')
        elif read.flag & 256:
            # use full junction sequence for assembly
            junctionfqfh.write('@' + read.query_name + '_' + str(read.flag) + '\n')
            junctionfqfh.write(read.query_sequence + '\n')
            junctionfqfh.write('+' + '\n')
            junctionfqfh.write(''.join(
                list(map(chr, [x + 33 for x in read.query_qualities]))) + '\n')
            # write overhang separately
            overhangfqfh.write('@' + read.query_name + '_' + str(read.flag) + '\n')
            overhangfqfh.write(read.query_alignment_sequence + '\n')
            overhangfqfh.write('+' + '\n')
            overhangfqfh.write(''.join(
                list(map(chr, [x + 33 for x in read.query_alignment_qualities]))) + '\n')
    pairfqfh.close()
    junctionfqfh.close()
    os.remove(bam_sort + '.bam')


def get_rna_support(jxn, tx, s_reads, in_bam, gtree):
    ''' Run command to identify read support for a single breakpoint '''
    logger.debug('Getting read support for ' + jxn)
    start = time.time()
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = 'support' + '/' + clean_jxn + '/'
    su.common.make_new_dir(jxn_dir)
    # get supporting read counts and other metrics
    s_reads = s_reads.split(',')
    extracted = get_reads_from_bam(in_bam, jxn, tx, s_reads, gtree)
    logger.debug('Found and extracted reads successfully')
    results = collections.defaultdict(list)
    results['name'] = jxn
    u_reads = []
    for ktype in extracted:
        for ktuple in extracted[ktype]:
            kstrand, korient = ktuple
            dname = ktype + '_' + kstrand + '_' + korient
            found_reads = list(set(extracted[ktype][ktuple].keys()))
            results[dname + "_reads"] = found_reads
            results[dname] = len(results[dname + "_reads"])
            for kname, val in six.iteritems(extracted[ktype][ktuple]):
                u_reads.append(kname)
                vun, vas, vmm, vseqlen, vmeanbq = val  # metric values: read flag, AS, nM, length, qualities
                results[dname + '_AS'].append(vas)
                results[dname + '_mismatches'].append(vmm)
                results[dname + '_seqlen'].append(vseqlen)
                results[dname + '_meanBQ'].append(vmeanbq)
            # convert lists to strings for easier record keeping
            results[dname + '_reads'] = ','.join(list(map(str, results[dname + '_reads'])))
            results[dname + '_AS'] = ','.join(list(map(str, results[dname + '_AS'])))
            results[dname + '_mismatches'] = ','.join(list(map(str, results[dname + '_mismatches'])))
            results[dname + '_seqlen'] = ','.join(list(map(str, results[dname + '_seqlen'])))
            results[dname + '_meanBQ'] = ','.join(list(map(str, results[dname + '_meanBQ'])))

    # write unique reads to file
    read_ids_all = jxn_dir + 'supporting_reads_all.txt'
    f2 = open(read_ids_all, 'w')
    for uread in set(u_reads):
        f2.write(uread + '\n')
    f2.close()

    # subset supporting reads to bam
    support_bam = jxn_dir + 'supporting.bam'
    logger.debug('Subsetting reads from bam')
    subset_bam_by_reads(in_bam, support_bam, read_ids_all, jxn)
    pysam.index(support_bam)
    # convert to fastq
    pairfq = jxn_dir + 'paired.fastq'
    junctionfq = jxn_dir + 'junctions.fastq'
    overhangfq = jxn_dir + 'overhang.fastq'
    bam2fastq(jxn_dir, support_bam, junctionfq, pairfq, overhangfq)
    logger.info(jxn + ' support took  %g seconds' % (time.time() - start))
    return results
