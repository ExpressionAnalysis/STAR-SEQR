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
import pysam  # requires 0.9.0 or newer
import multiprocessing as mp
import signal
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import starseqr_utils
import annotate_sv as ann

__author__ = "Jeff Jasper"
__email__ = "jasper1918@gmail.com"

logger = logging.getLogger("STAR-SEQR")


def make_new_dir(newdir):
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass


def find_resource(filename):
    packagedir = starseqr_utils.__path__[0]
    dirname = os.path.join(packagedir, 'resources')
    fullname = os.path.abspath(os.path.join(dirname, filename))
    return fullname


def find_discspan(jxn, bam):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pairDict = {}
    pairDict[('disc', 'all')] = {}
    pairDict[('disc', 'unique')] = {}
    # consider 129,65(F,F)  97,145(F,R) 113,177(R,R)
    # Can these ever be proper paired? What if short indel?
    # need the minus 1 to be exact.
    if str1 == "+":
        pos1left = int(pos1) - 100000
        pos1right = int(pos1) + int(repright) - 1
    if str2 == "-":
        pos2left = int(pos2) - 100000
        pos2right = int(pos2) + int(repright) - 1
    if str1 == "-":
        pos1left = int(pos1) - int(repright) - 1
        pos1right = int(pos1) + 100000
    if str2 == "+":
        pos2left = int(pos2) - int(repright) - 1
        pos2right = int(pos2) + 100000
    # print(pos1left, pos1right, pos2left, pos2right)
    for read in bam:
        if (read.next_reference_name == chrom2 and
                read.reference_name == chrom1 and
                not read.flag & 256):
            if (int(read.reference_start) > pos1left and
                    int(read.reference_start) < pos1right and
                    int(read.next_reference_start) > pos2left and
                    int(read.next_reference_start) < pos2right):
                if int(read.get_tag('AS')) > 1:  # filter on local alignment score
                    # print(read.reference_start, read.next_reference_start,
                    #  read.query_name, read.flag, read.mate_is_reverse)
                    if read.is_paired:
                        pairDict[('disc', 'all')][read.query_name] = read.flag
                    if not read.is_duplicate and read.is_paired:
                        pairDict[('disc', 'unique')][read.query_name] = read.flag
    return pairDict


def find_junctions(bam, chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright):
    '''
    flags: 321 and 385 are orient1. 337 and 401 are orient2.
    '''
    # chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    jxnDict = {}
    jxnDict[('all', 'for', 'first')] = {}
    jxnDict[('all', 'for', 'second')] = {}
    jxnDict[('all', 'rev', 'first')] = {}
    jxnDict[('all', 'rev', 'second')] = {}
    jxnDict[('unique', 'for', 'first')] = {}
    jxnDict[('unique', 'for', 'second')] = {}
    jxnDict[('unique', 'rev', 'first')] = {}
    jxnDict[('unique', 'rev', 'second')] = {}

    if str1 == "+":
        pos1left = int(pos1) - 10000
        pos1right = int(pos1) + int(repright) - 1  # -1 ?
    if str2 == "-":
        pos2left = int(pos2) - 10000
        pos2right = int(pos2) + int(repright) - 1
    if str1 == "-":
        pos1left = int(pos1) - int(repright) - 1
        pos1right = int(pos1) + 10000
    if str2 == "+":
        pos2left = int(pos2) - int(repright) - 1
        pos2right = int(pos2) + 10000
    # print(pos1left, pos1right, pos2left, pos2right)

    for read in bam:
        if (read.next_reference_name == chrom2 and
                read.reference_name == chrom1 and
                read.flag & 256):
            # print(read.reference_start, read.next_reference_start,
            #      read.query_name, read.flag, read.mate_is_reverse)
            if (int(read.reference_start) > pos1left and
                    int(read.reference_start) < pos1right or
                    int(read.next_reference_start) > pos2left and
                    int(read.next_reference_start) < pos2right):
                # print(read.query_name, read.flag, read.reference_name,
                #      read.reference_start, read.next_reference_name, read.next_reference_start)
                # todo: check that read falls on same gene as breakpoint.
                if int(read.get_tag('AS')) > 1:
                    if read.is_paired:
                        if read.flag & 16:
                            if read.flag & 64:
                                jxnDict[('all', 'rev', 'first')][read.query_name] = read.flag
                            else:
                                jxnDict[('all', 'rev', 'second')][read.query_name] = read.flag
                        else:
                            if read.flag & 64:
                                jxnDict[('all', 'for', 'first')][read.query_name] = read.flag
                            else:
                                jxnDict[('all', 'for', 'second')][read.query_name] = read.flag
                    if not read.is_duplicate and read.is_paired:
                        if read.flag & 16:
                            if read.flag & 64:
                                jxnDict[('unique', 'rev', 'first')][read.query_name] = read.flag
                            else:
                                jxnDict[('unique', 'rev', 'second')][read.query_name] = read.flag
                        else:
                            if read.flag & 64:
                                jxnDict[('unique', 'for', 'first')][read.query_name] = read.flag
                            else:
                                jxnDict[('unique', 'for', 'second')][read.query_name] = read.flag
    return jxnDict


def get_reads_from_bam(bam_file, jxn, args):
    '''
    Fetch reads from each direction of jxn and get an accounting of each if unique or duplicate.
    '''
    logger.debug("Extracting supporting read information for " + jxn)
    results = {}
    if not os.path.exists(bam_file):
        logger.error("BAM file could not be found: " + bam_file)
        sys.exit()
    bamObject = pysam.Samfile(bam_file, 'rb')
    if not bamObject.check_index():
        logger.error("Index file could not be found for: " + bam_file)
        sys.exit()
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    # account for genome start
    if int(pos1) - 1000 > 0:
        dist1 = int(pos1)
    else:
        dist1 = 0
    if int(pos2) - 1000 > 0:
        dist2 = int(pos2)
    else:
        dist2 = 0
    # account for genome boundaries
    # ++Todo: Need to get bam header to look for ends. write fxn for this. also remove circular chrom
    logger.debug("Extracting paired spanning reads for " + jxn)
    # spans
    bam = bamObject.fetch(chrom1, int(pos1) - int(dist1), int(pos1) + dist1)
    spanD = find_discspan(jxn, bam)
    results['spans'] = spanD
    # left junctions
    logger.debug("Extracting junction reads originating from the first breakpoint for " + jxn)
    bam = bamObject.fetch(chrom1, int(pos1) - int(dist1), int(pos1) + dist1)
    jxnD_for = find_junctions(bam, chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright)
    results['jxnleft'] = jxnD_for
    # right junctions
    logger.debug("Extracting junction reads originating from the second breakpoint for " + jxn)
    flipstr = string.maketrans("-+", "+-")
    nstr1 = str1.translate(flipstr)
    nstr2 = str2.translate(flipstr)
    bam = bamObject.fetch(chrom2, int(pos2) - int(dist2), int(pos2) + dist2)
    jxnD_rev = find_junctions(bam, chrom2, pos2, nstr2, chrom1, pos1, nstr1, repleft, repright)
    results['jxnright'] = jxnD_rev
    bamObject.close()
    logger.debug("Finished extracting supporting read information for " + jxn)
    return results


def subset_bam_by_reads(bam, out_bam, read_ids):
    logger.debug("Subset bam with supporting reads to " + out_bam)
    names = "names=" + read_ids
    index = "index=1"
    indexfilename = "indexfilename=" + out_bam + ".bai"
    tmpfile = "tmpfile=" + out_bam + ".tmp"
    stdinfile = bam
    bamf = open(out_bam, "wb")
    try:
        with open(stdinfile, "rb") as f:
            retcode = sp.call(['bamfilternames', names, index, indexfilename,
                               tmpfile], stdin=f, stdout=bamf, stderr=sp.PIPE)
        bamf.close()
        if retcode != 0:
            logger.error("bamfilternames failed!")
            sys.exit(1)
    except OSError, o:
        logger.error("Exception: " + str(o))
        logger.error("bamfilternames Failed", exc_info=True)
        bamf.close()
        sys.exit(1)


def bam_2_nsort_bam(in_bam):
    bam_sort = in_bam[:-4] + ".nsorted"
    pysam.sort("-n", in_bam, "-o", bam_sort + ".bam")


# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}


def reverse_complement(sequence):
    return ''.join(COMPLEMENT[x] for x in sequence[::-1])


def bam2fastq(jxn_dir, in_bam, junctionfq, pairfq):
    '''
    STAR writes 3 lines in the SAM file for chimeric reads.
    Custom convert bam2fastq where pair1/2 go into paired.fastq
    and chimeric only reads go into junctions.fq.
    Chimeric reads need all seqs/quals while
    paired must be soft-clipped seqs/quals.
    '''
    logger.debug("Converting bam to fastq for " + in_bam)
    if not os.path.exists(in_bam):
        logger.error("BAM file could not be found: " + in_bam)
        sys.exit(1)
    pairfqfh = open(pairfq, 'w')
    junctionfqfh = open(junctionfq, 'w')
    bam_sort = in_bam[:-4] + ".nsorted"
    pysam.sort("-n", in_bam, "-o", bam_sort + ".bam")
    bamObject = pysam.Samfile(bam_sort + ".bam", 'rb')
    for read in bamObject.fetch(until_eof=True):
        if not read.flag & 256:
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
            pairfqfh.write("@" + read.query_name + "/" + str(orient) + "\n")
            pairfqfh.write(seq + "\n")
            pairfqfh.write("+" + "\n")
            pairfqfh.write(''.join(
                map(chr, [x + 33 for x in quals])) + "\n")
        elif read.flag & 256:
            junctionfqfh.write("@" + read.query_name + "_" + str(read.flag) + "\n")
            junctionfqfh.write(read.query_sequence + "\n")
            junctionfqfh.write("+" + "\n")
            junctionfqfh.write(''.join(
                map(chr, [x + 33 for x in read.query_alignment_qualities])) + "\n")
    os.remove(bam_sort + ".bam")


def fasta_iter(fasta_name):
    """
    Given a fasta file, yield tuples of header, sequence
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        # print(seq)
        yield header, seq


def do_velvet(assemdir, fastq, kmer, errlog, *args):
    '''
    Run velvet with single or multiple fastqs.
    '''
    logger.debug("*Running Velvet.")
    velveth_cmd = ["velveth", assemdir, str(kmer), "-short", "-fastq", fastq]
    for ar in args:
        velveth_cmd.extend(["-shortPaired", "-fastq", ar])
    vh_args = map(str, velveth_cmd)
    logger.debug("*velveth Command: " + " ".join(vh_args))
    try:
        p = sp.Popen(velveth_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error("Error: velveth failed")
            sys.exit(1)
    except (OSError) as o:
        logger.error("Exception: " + str(o))
        logger.error("velveth Failed!", exc_info=True)
        sys.exit(1)
    velvetg_cmd = ["velvetg", assemdir, "-cov_cutoff", "2"]
    vg_args = map(str, velvetg_cmd)
    logger.debug("*velvetg Command: " + " ".join(vg_args))
    try:
        p = sp.Popen(velvetg_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error("Error: velvetg failed")
            sys.exit(1)
    except (OSError) as o:
        logger.error("Exception: " + str(o))
        logger.error("velvetg Failed!", exc_info=True)
        sys.exit(1)
    # extract header, sequence to tuple
    records = fasta_iter(assemdir + "/contigs.fa")
    # for x in records:
    #     print(x)
    # just return first record for now
    # return islice(records, 1)
    return list(islice(records, 0, 1))


def do_spades(assemdir, pfastq, jxnfastq, errlog):
    logger.debug("*Running SPADES")
    spades_cmd = ['spades.py', "--12", pfastq, "-s", jxnfastq,
                  "-o", assemdir, "--careful", "--phred-offset", 33,
                  "-t", "1", "-m", "5",
                  "-cov-cutoff", "off"]
    spades_args = map(str, spades_cmd)
    logger.debug("*SPADES Command: " + " ".join(spades_args))
    try:
        p = sp.Popen(spades_args, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            errlog.write(stdout)
        if stderr:
            errlog.write(stderr)
        if p.returncode != 0:
            logger.error("Error: spades failed")
            # sys.exit(1)
        else:
            if (os.stat(os.path.realpath(assemdir + "/scaffolds.fasta")).st_size != 0):
                records2 = fasta_iter(assemdir + "/scaffolds.fasta")
                return list(islice(records2, 0, 1))
    except (OSError) as o:
        logger.error("Exception: " + str(o))
        logger.error("SPADES Failed!", exc_info=True)
        sys.exit(1)


def run_support_fxn(jxn, in_bam, args, *opts):
    ''' Run command to run assembly for a single breakpoint '''
    logger.info("Getting read support for " + jxn)
    start = time.time()
    results = {}
    clean_jxn = str(jxn).replace(":", "_")
    clean_jxn = str(clean_jxn).replace("+", "pos")
    clean_jxn = str(clean_jxn).replace("-", "neg")
    jxn_dir = "support" + "/" + clean_jxn + "/"
    make_new_dir(jxn_dir)
    # get supporting reads
    extracted = get_reads_from_bam(in_bam, jxn, args)
    logger.debug("Found and extracted reads successfully")
    for entry in extracted:
        for entry2 in extracted[entry]:
            # print(entry + '_' + '_'.join(entry2), len(extracted[entry][entry2].keys()))
            results[entry + '_' + '_'.join(entry2)] = len(extracted[entry][entry2].keys())
    # unique
    read_ids = jxn_dir + "supporting_reads_unique.txt"
    with open(read_ids, "w") as f:
        for entry in extracted:
            for entry2 in extracted[entry]:
                if "unique" in entry2:
                    for uniq_reads in extracted[entry][entry2]:
                        f.write(uniq_reads + "\n")
    f.close
    # all_reads
    read_ids_all = jxn_dir + "supporting_reads_all.txt"
    with open(read_ids_all, "w") as f2:
        for entry in extracted:
            for entry2 in extracted[entry]:
                for all_reads in extracted[entry][entry2]:
                    f2.write(all_reads + "\n")
    f2.close
    # subset supporting reads to bam
    support_bam = jxn_dir + "supporting.bam"
    logger.debug("*Subsetting reads from bam")
    subset_bam_by_reads(in_bam, support_bam, read_ids)
    pysam.index(support_bam)
    # convert to fastq and run velvet
    pairfq = jxn_dir + 'paired.fastq'
    junctionfq = jxn_dir + 'junctions.fastq'
    bam2fastq(jxn_dir, support_bam, junctionfq, pairfq)
    errlog = open(jxn_dir + "assembly_log.txt", "w")
    # velvet_all = do_velvet(jxn_dir + "assem_pair", junctionfq, 17, errlog, pairfq)
    velvet_all = do_velvet(jxn_dir + "assem_jxn", junctionfq, 17, errlog)
    errlog.close()
    if velvet_all:
        vname, vseq = velvet_all[0]
        results['velvet'] = vseq
    if args.spades:
        splog = open(jxn_dir + "spades_log.txt", "w")
        spades_seq = do_spades(jxn_dir + "spades", pairfq, junctionfq, splog)
        splog.close()
        if spades_seq:
            sname, sseq = spades_seq[0]
            results['spades'] = sseq
    results['name'] = jxn
    end = time.time()
    elapsed = end - start
    logger.info("Support took  %g seconds" % (elapsed))
    return results


def init_worker():
    '''this is to sidestep multiproc bug with keyboard interrupts'''
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def run_support_parallel(df, in_bam, args):
    '''Run assembly in parallel'''
    logger.debug("Getting Read Support")
    make_new_dir("support")
    results = []
    pool = mp.Pool(int(args.workers), init_worker)
    try:
        for jxn in df['name']:
            seq = pool.apply_async(run_support_fxn, args=[jxn, in_bam, args])
            results.append(seq)
            # results[jxn] = jxnres
        pool.close()
    except KeyboardInterrupt as e:
        logger.error("Error: Keyboard interrupt")
        pool.terminate()
        raise e
    except Exception as e:
        pool.terminate()
        raise e
    finally:
        pool.join()
    list_of_dicts = []
    for res in results:
        jxn_ld = res.get()
        list_of_dicts.append(jxn_ld)
        # for x in jxn_ld:
        #     print(x,jxn_ld[x])
    logger.debug("Finished Assembly")
    return list_of_dicts
