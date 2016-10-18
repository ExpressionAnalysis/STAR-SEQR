#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import re
import string
import time
import gzip
import errno
import logging
import subprocess as sp
from itertools import groupby, islice
import pysam  # requires 0.9.0
import multiprocessing as mp
import signal
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import collections
import numpy as np


logger = logging.getLogger('STAR-SEQR')


def make_new_dir(newdir):
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass


def find_resource(filename):
    packagedir = os.path.dirname(os.path.realpath(__file__))
    dirname = os.path.join(packagedir, 'resources')
    fullname = os.path.abspath(os.path.join(dirname, filename))
    return fullname


def get_pos_genes(chrom1, pos1, gtree):
    resL = gtree[chrom1].search(int(pos1))
    genesL = set()
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            Lsymbol = list(resL)[idx].data['name2']
            genesL.add(Lsymbol)
    else:
        genesL.add("NA")
    return list(genesL)


def find_support_reads(jxn, bam, side, tx, sub_reads, gtree):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1 = int(pos1)
    pos2 = int(pos2)

    if side == 2:
        chrom1, chrom2 = chrom2, chrom1
        pos1, pos2 = pos2, pos1
        flipstr = string.maketrans('-+', '+-')
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
                read_tx = get_pos_genes(read.reference_name, read.reference_start, gtree)
                read_tx2 = get_pos_genes(read.next_reference_name, read.next_reference_start, gtree)
                if set(tx).intersection(set(read_tx)) and set(tx).intersection(set(read_tx2)):
                    retDict[mytype][(mystrand, myorient)][read.query_name] = (is_subset,
                                                                              read.get_tag('AS'),
                                                                              read.get_tag('nM'),
                                                                              read.query_alignment_length,
                                                                              int(np.mean(read.query_alignment_qualities)))
    return retDict


def get_reads_from_bam(bam_file, jxn, tx, s_reads, gtree, args):
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
    len1 = (item for item in bamObject.header['SQ'] if item['SN'] == chrom1).next()['LN']
    len2 = (item for item in bamObject.header['SQ'] if item['SN'] == chrom2).next()['LN']
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
            # !!! Don't do a sys.exit here.. you will hang and lose error info
            logger.error('bamfilternames failed on:' + jxn)
            # for line in retcode.stderr:
            #     logger.error(line)
            # os._exit(1)
    except OSError, o:
        logger.error('bamfilternames Failed', exc_info=True)
        logger.error('Exception: ' + str(o))
        bamf.close()
        os._exit(1)


def bam_2_nsort_bam(in_bam):
    bam_sort = in_bam[:-4] + '.nsorted'
    pysam.sort('-n', in_bam, '-o', bam_sort + '.bam')


def reverse_complement(sequence):
    # DNA base complements
    COMPLEMENT = {'A': 'T',
                  'T': 'A',
                  'C': 'G',
                  'G': 'C',
                  'N': 'N'}
    return ''.join(COMPLEMENT[x] for x in sequence[::-1])


def bam2fastq(jxn_dir, in_bam, junctionfq, pairfq):
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
                    seq = reverse_complement(read.query_alignment_sequence)
                    quals = read.query_alignment_qualities[::-1]
                else:
                    seq = read.query_alignment_sequence
                    quals = read.query_alignment_qualities
            else:
                orient = 2
                if read.is_reverse:
                    seq = reverse_complement(read.query_alignment_sequence)
                    quals = read.query_alignment_qualities[::-1]
                else:
                    seq = read.query_alignment_sequence
                    quals = read.query_alignment_qualities
            pairfqfh.write('@' + read.query_name + '/' + str(orient) + '\n')
            pairfqfh.write(seq + '\n')
            pairfqfh.write('+' + '\n')
            pairfqfh.write(''.join(
                map(chr, [x + 33 for x in quals])) + '\n')
        elif read.flag & 256:
            junctionfqfh.write('@' + read.query_name + '_' + str(read.flag) + '\n')
            junctionfqfh.write(read.query_sequence + '\n')
            junctionfqfh.write('+' + '\n')
            junctionfqfh.write(''.join(
                map(chr, [x + 33 for x in read.query_alignment_qualities])) + '\n')
    pairfqfh.close()
    junctionfqfh.close()
    os.remove(bam_sort + '.bam')


def fasta_iter(fasta_name):
    '''
    Given a fasta file, yield tuples of header, sequence
    '''
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
    '''
    Run velvet with single or multiple fastqs.
    '''
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
            sys.exit(1)
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
    # extract header, sequence to tuple
    records = fasta_iter(assemdir + '/contigs.fa')
    # for x in records:
    #     print(x)
    # just return first record for now
    # return islice(records, 1)
    return list(islice(records, 0, 1))


def do_spades(assemdir, pfastq, jxnfastq, errlog):
    logger.debug('*Running SPADES')
    spades_cmd = ['spades.py', '--12', pfastq, '-s', jxnfastq,
                  '-o', assemdir, '--careful', '--phred-offset', 33,
                  '-t', '1', '-m', '5',
                  '-cov-cutoff', 'off']
    spades_args = map(str, spades_cmd)
    logger.debug('*SPADES Command: ' + ' '.join(spades_args))
    try:
        p = sp.Popen(spades_args, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error('Error: spades failed')
            # sys.exit(1)
        else:
            if (os.stat(os.path.realpath(assemdir + '/scaffolds.fasta')).st_size != 0):
                records2 = fasta_iter(assemdir + '/scaffolds.fasta')
                return list(islice(records2, 0, 1))
    except (OSError) as o:
        logger.error('Exception: ' + str(o))
        logger.error('SPADES Failed!', exc_info=True)
        sys.exit(1)


def run_support_fxn(jxn, tx, s_reads, in_bam, args, *opts):
    ''' Run command to run assembly for a single breakpoint '''
    logger.debug('Getting read support for ' + jxn)
    start = time.time()
    clean_jxn = str(jxn).replace(':', '_')
    clean_jxn = str(clean_jxn).replace('+', 'pos')
    clean_jxn = str(clean_jxn).replace('-', 'neg')
    jxn_dir = 'support' + '/' + clean_jxn + '/'
    make_new_dir(jxn_dir)
    # get supporting read counts and other metrics
    extracted = get_reads_from_bam(in_bam, jxn, tx, s_reads, gtree, args)
    logger.debug('Found and extracted reads successfully')

    results = collections.defaultdict(list)
    u_reads = []
    for ktype in extracted:
        for ktuple in extracted[ktype]:
            kstrand, korient = ktuple
            dname = ktype + '_' + kstrand + '_' + korient
            found_reads = list(set(extracted[ktype][ktuple].keys()))
            results[dname + "_reads"] = found_reads
            results[dname] = len(results[dname + "_reads"])
            for kname, val in extracted[ktype][ktuple].iteritems():
                u_reads.append(kname)
                vun, vas, vmm, vseqlen, vmeanbq = val  # metric values: read flag, AS, nM, length, qualities
                results[dname + '_AS'].append(vas)
                results[dname + '_mismatches'].append(vmm)
                results[dname + '_seqlen'].append(vseqlen)
                results[dname + '_meanBQ'].append(vmeanbq)

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
    # convert to fastq and run velvet
    pairfq = jxn_dir + 'paired.fastq'
    junctionfq = jxn_dir + 'junctions.fastq'
    bam2fastq(jxn_dir, support_bam, junctionfq, pairfq)
    errlog = open(jxn_dir + 'assembly_log.txt', 'w')
    velvet_all = do_velvet(jxn_dir + 'assem_pair', junctionfq, 17, errlog, pairfq)
    # velvet_all = do_velvet(jxn_dir + 'assem_jxn', junctionfq, 17, errlog)
    errlog.close()
    if velvet_all:
        vname, vseq = velvet_all[0]
        results['assembly'] = vseq
    elif args.spades:
        splog = open(jxn_dir + 'spades_log.txt', 'w')
        spades_seq = do_spades(jxn_dir + 'spades', pairfq, junctionfq, splog)
        splog.close()
        if spades_seq:
            sname, sseq = spades_seq[0]
            results['assembly'] = sseq
    else:
        results['assembly'] = ''
    results['name'] = jxn
    logger.info(jxn + ' support took  %g seconds' % (time.time() - start))
    return results


def init_worker():
    '''this is to sidestep multiproc bug with keyboard interrupts'''
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def run_support_parallel(df, in_bam, args):
    '''Run assembly in parallel'''
    logger.debug('Getting Read Support')
    make_new_dir('support')
    results = []
    # annotation to be used
    global gtree
    if args.ann_source == "refgene":
        refgene = find_resource("refGene.txt.gz")
        kg = gzip.open(refgene)
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.REF_GENE)
    elif args.ann_source == "ensgene":
        ensgene = find_resource("ensGene.txt.gz")
        kg = gzip.open(ensgene)
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)
    elif args.ann_source == "gencode":
        gencode = find_resource("wgEncodeGencodeBasicV24lift37.txt.gz")
        kg = gzip.open(gencode)
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)

    pool = mp.Pool(int(args.threads), init_worker)
    try:
        for i in (df.index.values):
            args2pass = [df.loc[i, 'name'], df.loc[i, 'txunion'], df.loc[i, 'supporting_reads'].split(','), in_bam, args]
            seq = pool.apply_async(run_support_fxn, args=args2pass)
            results.append(seq)
        pool.close()
        pool.join()
        list_of_dicts = []
        for res in results:
            jxn_ld = res.get()
            list_of_dicts.append(jxn_ld)
    except KeyboardInterrupt as e:
        logger.error('Error: Keyboard interrupt')
        pool.terminate()
        raise e
    except Exception as e:
        logger.error('Exception: ' + str(e))
        pool.terminate()
        raise e

    logger.debug('Finished Assembly')
    return list_of_dicts
