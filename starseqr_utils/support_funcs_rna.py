#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import six
from six import reraise as raise_
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


def find_support_reads(jxn, bam, side, tx, sub_reads, gtree, chimflag):
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
                not (read.flag & 2) and not (read.flag & chimflag)):
            if ((int(read.reference_start) > pos1left and
                 int(read.reference_start) < pos1right) and
                (int(read.next_reference_start) > pos2left and
                 int(read.next_reference_start) < pos2right)):
                mymode = 'span'
        elif (read.next_reference_name == read.reference_name and
                (read.flag & 2) and not (read.flag & chimflag)):
            if ((int(read.reference_start) > pos1left and
                 int(read.reference_start) < pos1right) or
                (int(read.next_reference_start) > pos2left and
                 int(read.next_reference_start) < pos2right)):
                mymode = 'jxn'
        elif (read.next_reference_name == chrom2 and read.reference_name == chrom1 and
                not (read.flag & 2) and (read.flag & chimflag)):
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


def get_reads_from_bam(bam_file, jxn, tx, s_reads, gtree, chimflag):
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
    support_for = find_support_reads(jxn, bam, 1, tx, s_reads, gtree, chimflag)
    results['spanleft'] = support_for['spanleft']
    results['jxnleft'] = support_for['jxnleft']
    results['hangleft'] = support_for['hangleft']
    # support right
    logger.debug('Extracting right paired reads for ' + jxn)
    bam = bamObject.fetch(chrom2, dist2_less, dist2_plus)
    support_rev = find_support_reads(jxn, bam, 2, tx, s_reads, gtree, chimflag)
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
            sp.call(['bamfilternames', names, index, indexfilename, tmpfile], stdin=f, stdout=bamf, stderr=sp.PIPE)
        bamf.close()
    except OSError as e:
        bamf.close()
        logger.error('bamfilternames Failed', exc_info=True)
        logger.error('Exception: ' + str(e))
        traceback = sys.exc_info()[2]
        raise_(ValueError, e, traceback)
    return


def bam2fastq(in_bam, bam_type, jxn_dir, chimflag):
    logger.debug('Converting bams to fastqs')
    if not os.path.exists(in_bam):
        logger.error('BAM file could not be found: ' + in_bam)
        sys.exit(1)
    try:
        bam_nsort = in_bam[:-4] + '.nsorted.bam'
        pysam.sort('-n', in_bam, '-o', bam_nsort, catch_stdout=False)
    except (OSError) as e:
        logger.error('Exception: ' + str(e))
        logger.error('pysam Failed!', exc_info=True)
        traceback = sys.exc_info()[2]
        raise_(ValueError, e, traceback)
    bamObject = pysam.Samfile(bam_nsort, 'rb')

    if bam_type == "span":
        out_fq = open(jxn_dir + 'span.fastq', 'w')
        for read in bamObject.fetch(until_eof=True):
            # use aligned sequenc for spanning reads
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
            out_fq.write('@' + read.query_name + '/' + str(orient) + '\n')
            out_fq.write(seq + '\n')
            out_fq.write('+' + '\n')
            out_fq.write(''.join(
                list(map(chr, [x + 33 for x in quals]))) + '\n')
        out_fq.close()
        os.remove(bam_nsort)
    elif bam_type == "split":
        out_fq = open(jxn_dir + 'split.fastq', 'w')
        out_fq2 = open(jxn_dir + 'overhang.fastq', 'w')
        for read in bamObject.fetch(until_eof=True):
            # use full junction sequence for split reads...accomodates assembly
            if not read.flag & chimflag:
                out_fq.write('@' + read.query_name + '_' + str(read.flag) + '\n')
                out_fq.write(read.query_sequence + '\n')
                out_fq.write('+' + '\n')
                out_fq.write(''.join(
                    list(map(chr, [x + 33 for x in read.query_qualities]))) + '\n')
            # write overhang separately
            elif read.flag & chimflag:
                out_fq2.write('@' + read.query_name + '_' + str(read.flag) + '\n')
                out_fq2.write(read.query_alignment_sequence + '\n')
                out_fq2.write('+' + '\n')
                out_fq2.write(''.join(
                    list(map(chr, [x + 33 for x in read.query_alignment_qualities]))) + '\n')
        out_fq.close()
        out_fq2.close()
        os.remove(bam_nsort)
    return


def starbam2genesupport(in_bam_path, norm_bam_path, chim_bam_path, region1, region2):
    # get chimeric reads using ch tag from "WithinBAM method", split into categories
    # You will have errors if resulting bams are not sorted!
    infile = pysam.AlignmentFile(in_bam_path, "rb")
    # norm_out = pysam.AlignmentFile(norm_bam_path + 'unsorted.bam', "wb", template=infile)
    chim_out = pysam.AlignmentFile(chim_bam_path + 'unsorted.bam', "wb", template=infile)
    for r in infile.fetch(region=region1):
        if r.has_tag("ch"):
            chim_out.write(r)
        # if not r.flag & 2048:
        #     norm_out.write(r)
    for r in infile.fetch(region=region2):
        if r.has_tag("ch"):
            chim_out.write(r)
        # if not r.flag & 2048:
        #     norm_out.write(r)
    infile.close()
    # norm_out.close()
    chim_out.close()
    # su.common.sortnindex_bam(norm_bam_path + 'unsorted.bam', norm_bam_path, 1)
    su.common.sortnindex_bam(chim_bam_path + 'unsorted.bam', chim_bam_path, 1)
    return


def normbamtofastq(bam, jxn_dir):
    # bamfilternames is sensitive to duplicate read names in the list. Make sure they are unique.
    logger.debug('Subset bam with supporting reads to fastqs for salmon')
    filename = 'filename=' + bam
    f1 = 'F=' + jxn_dir + 'read1.fq'
    f2 = 'F2=' + jxn_dir + 'read2.fq'
    o1 = 'O=' + jxn_dir + 'un1.fq'
    o2 = 'O2=' + jxn_dir + 'un2.fq'
    cmdargs = ['bamtofastq', filename, 'inputformat=bam', f1, f2, o1, o2]
    # logger.info("*Command: " + " ".join(cmdargs))
    try:
        p = sp.Popen(cmdargs, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            pass
            # logger.info(stdout)
        if stderr:
            pass
            # logger.error(stderr)
    except OSError as e:
        logger.error('bamtofastq Failed', exc_info=True)
        logger.error('Exception: ' + str(e))
        traceback = sys.exc_info()[2]
        raise_(ValueError, e, traceback)
    return


def get_rna_support(jxn, tx, s_reads, in_bam, gtree, chimflag):
    ''' Run command to identify read support for a single breakpoint '''
    logger.debug('Getting read support for ' + jxn)
    start = time.time()
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = 'support' + '/' + clean_jxn + '/'
    su.common.make_new_dir(jxn_dir)
    # initialize results
    results = collections.defaultdict(list)
    results['name'] = jxn
    # get genes bounds and extract all reads to subset bams
    if chimflag == 2048:
        region1 = su.annotate_sv.get_gene_region(chrom1, pos1, gtree)
        region2 = su.annotate_sv.get_gene_region(chrom2, pos2, gtree)
        normbam = jxn_dir + 'normal.bam'
        chimbam = jxn_dir + 'chimeric.bam'
        starbam2genesupport(in_bam, normbam, chimbam, region1, region2)
        # normbamtofastq(normbam, jxn_dir)
        in_bam = chimbam
    # get supporting read counts and other metrics
    s_reads = s_reads.split(',')
    extracted = get_reads_from_bam(in_bam, jxn, tx, s_reads, gtree, chimflag)
    logger.debug('Found and extracted reads successfully')
    span_reads = []
    split_reads = []
    for ktype in extracted:
        for ktuple in extracted[ktype]:
            kstrand, korient = ktuple
            dname = ktype + '_' + kstrand + '_' + korient
            found_reads = list(set(extracted[ktype][ktuple].keys()))
            results[dname + "_reads"] = found_reads
            results[dname] = len(results[dname + "_reads"])
            for kname, val in six.iteritems(extracted[ktype][ktuple]):
                if str(ktype).startswith("span"):
                    span_reads.append(kname)
                elif str(ktype).startswith("jxn"):
                    split_reads.append(kname)
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

    # write unique reads to file for subsetting bam later
    # TODO: Consider making one bam and dividing after.
    # span reads
    spanread_path = jxn_dir + 'span_readids.txt'
    span_fh = open(spanread_path, 'w')
    for uread in set(span_reads):
        span_fh.write(uread + '\n')
    span_fh.close()
    # split reads
    splitread_path = jxn_dir + 'split_readids.txt'
    split_fh = open(splitread_path, 'w')
    for uread in set(split_reads):
        split_fh.write(uread + '\n')
    split_fh.close()

    # subset supporting reads to bam
    logger.debug('Subsetting reads from bam')
    # span
    spanbam = jxn_dir + 'span.bam'
    subset_bam_by_reads(in_bam, spanbam, spanread_path, jxn)
    pysam.index(spanbam)
    # split
    splitbam = jxn_dir + 'split.bam'
    subset_bam_by_reads(in_bam, splitbam, splitread_path, jxn)
    pysam.index(splitbam)

    # convert bam support to fastqs
    bam2fastq(spanbam, "span", jxn_dir, chimflag)
    bam2fastq(splitbam, "split", jxn_dir, chimflag)  # also produces overhang.fastq

    logger.info(jxn + ' support took  %g seconds' % (time.time() - start))
    return results
