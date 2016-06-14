#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import re
import string
import errno
import logging
import subprocess as sp
from itertools import groupby, islice
import pysam  # requires 0.9.0 or newer
import multiprocessing as mp
import signal

__author__ = "Jeff Jasper"
__email__ = "jasper1918@gmail.com"

logger = logging.getLogger("STAR-SEQR")
print(pysam.__version__)


def check_file_exists(path):
    if (os.stat(os.path.realpath(path)).st_size == 0):
        logger.error("Exiting. Cannot find file: " + os.path.realpath(path))
        sys.exit(1)


def run_star(cfg, fq1, fq2, args):
    '''
    Run STAR alignment for RNA or DNA with different sensitivity parameters.
    '''
    logger.info("Starting STAR Alignment")
    if not os.path.isfile(args.prefix + ".Chimeric.out.junction"):
        if args.nucleic_type == "DNA":
            STAR_args = [cfg['star_exec'], '--readFilesIn', fq1, fq2, '--readFilesCommand', 'zcat',
                         '--runThreadN', args.threads, '--genomeDir', cfg['star_index'],
                         '--outFileNamePrefix ', args.prefix + ".", '--outSAMtype', 'None',
                         '--alignIntronMax', 1, '--alignMatesGapMax', 0, '--sjdbOverhang', 0,
                         '--chimOutType', 'SeparateSAMold', '--chimScoreJunctionNonGTAG', 0]
            # choose sensitivity mode
            if (args.mode == 0):
                sens_params = ['--chimSegmentMin', 10, '--chimJunctionOverhangMin', 10,
                               '--chimScoreMin', 1, '--chimScoreDropMax', 20,
                               '--chimScoreSeparation', 10, '--chimSegmentReadGapMax', 0,
                               '--chimFilter', 'None']  # , "--twopassMode", "None"]
            elif (args.mode == 1):
                sens_params = ['--chimSegmentMin', 5, '--chimJunctionOverhangMin', 5,
                               '--chimScoreMin', 0, '--chimScoreDropMax', 10,
                               '--chimScoreSeparation', 5, '--chimSegmentReadGapMax', 0,
                               '--chimFilter', 'None']
                # no affect
                # , '--outSJfilterCountTotalMin', '1 1 1 1',
                # '--outSJfilterCountUniqueMin', '1 1 1 1', '--outSJfilterOverhangMin', '12 12 12 12',
                # '--outSJfilterIntronMaxVsReadN', '10000 20000 30000', '--alignSJstitchMismatchNmax', '0 0 0 0']
            STAR_args.extend(sens_params)
            # Need to convert all to string
            STAR_args = map(str, STAR_args)
        elif args.nucleic_type == "RNA":
            STAR_args = [cfg['star_exec'], '--readFilesIn', fq1, fq2, '--readFilesCommand', 'zcat',
                         '--runThreadN', args.threads, '--genomeDir', cfg['star_index_rna'],
                         '--outFileNamePrefix ', args.prefix + ".", '--outSAMtype', 'None',
                         '--alignIntronMax', 200000, '--alignMatesGapMax', 200000,
                         '--chimOutType', 'SeparateSAMold', '--chimScoreJunctionNonGTAG', -1,
                         '--alignSJDBoverhangMin', 10]  # '--alignSJstitchMismatchNmax', '5 -1 5 5'] throws errors
            # choose sensitivity mode
            if (args.mode == 0):
                sens_params = ['--chimSegmentMin', 10, '--chimJunctionOverhangMin', 10,
                               '--chimScoreMin', 1, '--chimScoreDropMax', 20,
                               '--chimScoreSeparation', 10, '--chimSegmentReadGapMax', 3,
                               '--chimFilter', 'None']
            elif (args.mode == 1):
                sens_params = ['--chimSegmentMin', 5, '--chimJunctionOverhangMin', 5,
                               '--chimScoreMin', 0, '--chimScoreDropMax', 10,
                               '--chimScoreSeparation', 5, '--chimSegmentReadGapMax', 3,
                               '--chimFilter', 'None']
            STAR_args.extend(sens_params)
            # Need to convert all to string
            STAR_args = map(str, STAR_args)
        else:
            logger.error("Need to define nucleic_type as either RNA or DNA")
            sys.exit(1)
        logger.info("*STAR Command: " + " ".join(STAR_args))
        # run STAR
        try:
            p = sp.Popen(STAR_args, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()
            if stdout:
                logger.info(stdout)
            if stderr:
                logger.error(stderr)
            if p.returncode != 0:
                logger.error("Error: STAR failed", exc_info=True)
                sys.exit(1)
        except (OSError) as o:
            logger.error("Exception: " + str(o))
            logger.error("STAR Failed", exc_info=True)
            sys.exit(1)
        check_file_exists(args.prefix + ".Chimeric.out.junction")
        check_file_exists(args.prefix + ".Chimeric.out.sam")
    else:
        check_file_exists(args.prefix + ".Chimeric.out.junction")
        check_file_exists(args.prefix + ".Chimeric.out.sam")
        logger.warn("Skipping STAR alignment as files already exist!")
    logger.info("STAR Alignment Finished!")


def sam_2_coord_bam(in_sam, out_bam):
    if (out_bam[-4:] == ".bam"):
        bam_prefix = out_bam[:-4]
    else:
        bam_prefix = out_bam
    bam_unsort = bam_prefix + ".unsorted.bam"
    pysam.view("-Sbu", "-o%s" % bam_unsort, in_sam)
    pysam.sort(bam_unsort, "-o", bam_prefix + ".bam")
    pysam.index("%s.bam" % bam_prefix)
    os.remove(bam_unsort)


def bam_2_nsort_sam(in_bam, out_sam):
    bam_sort = in_bam[:-4] + ".nsorted"
    pysam.sort("-n", in_bam, "-o", bam_sort + ".bam")
    pysam.view("-h", "-o%s" % out_sam, (bam_sort + ".bam"))
    os.remove(bam_sort + ".bam")


def run_biobambam2_rmdup(in_bam, out_bam, cfg):
    biobambam2 = (cfg['biobambam2_bin'] + '/bammarkduplicates2')
    mark_samdups_cmd = [biobambam2, "=".join(["I", in_bam]), "=".join(["O", out_bam]), "=".join(["markthreads",
                        cfg['biobambam_threads']]), "=".join(["M", "Dup_metrics.txt"]), "=".join(["level", "1"]),
                        "=".join(["verbose", "0"])]
    rmdup_args = map(str, mark_samdups_cmd)
    logger.info("MarkDups Command: " + " ".join(rmdup_args))
    try:
        p = sp.Popen(mark_samdups_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            logger.info(stdout)
        if stderr:
            logger.debug(stderr)
        if p.returncode != 0:
            logger.error("Error: bammarkduplicates failed", exc_info=True)
            sys.exit(1)
    except OSError, o:
        logger.error("Exception: " + str(o))
        logger.error("Biobambam2 Failed!", exc_info=True)
        sys.exit(1)


def fix_chimeric_flags(in_sam, out_sam):
    '''
    Biobambam2 does not mark the chimeric fragment. STAR has 3 ids for chimeric reads.
    Need to hold all reads with same id to mark the fragment.
    '''
    fixSAM = open(out_sam, 'w')
    data = open(in_sam, 'r')
    prevName = ""
    readDict = {}
    dup = 0
    for read in data:
        if read[0] == '@':
            print(read.strip(), file=fixSAM)
            continue
        readList = read.strip().split('\t')
        # print jxnName
        readName = readList[0]
        if (readName != prevName and prevName != ""):  # hits a new block of ids
            for flag in readDict.keys():
                if int(flag) > 1024:
                    dup = 1
            if (dup == 1):
                for iread in readDict:
                    rList = readDict[iread]
                    if int(rList[1]) < 1024:
                        rList[1] = int(rList[1]) + 1024
                    print(*rList, sep="\t", file=fixSAM)
            else:
                for iread in readDict:
                    rList = readDict[iread]
                    print(*rList, sep="\t", file=fixSAM)
            # reset parameters, add new read
            dup = 0
            readDict = {}
            readDict[readList[1]] = readList
        else:
            readDict[readList[1]] = readList
        prevName = readName
    # write last block
    for flag in readDict.keys():
        if int(flag) > 1024:
            dup = 1
        if (dup == 1):
            for iread in readDict:
                rList = readDict[iread]
                if int(rList[1]) < 1024:
                    rList[1] = int(rList[1]) + 1024
                print(*rList, sep="\t", file=fixSAM)
        else:
            for iread in readDict:
                rList = readDict[iread]
                print(*rList, sep="\t", file=fixSAM)


def markdups(in_sam, out_bam, cfg):
    # get biobambam2 to cfg or param
    logger.info("Marking duplicate reads")
    sam_2_coord_bam(in_sam, "primary.bam")
    run_biobambam2_rmdup("primary.bam", "mrkdup_tmp.bam", cfg)
    check_file_exists("mrkdup_tmp.bam")
    os.remove("primary.bam")
    os.remove("primary.bam.bai")
    bam_2_nsort_sam("mrkdup_tmp.bam", "mrkdup_tmp.sam")
    os.remove("mrkdup_tmp.bam")
    fix_chimeric_flags("mrkdup_tmp.sam", "final_tmp.sam")
    os.remove("mrkdup_tmp.sam")
    sam_2_coord_bam("final_tmp.sam", out_bam)
    os.remove("final_tmp.sam")
    check_file_exists(out_bam)
    logger.info("Finished marking duplicate reads")


def get_genome_bounds(bam):
    pass


def find_discspan(jxn, bam):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pairDict = {}
    pairDict[('all', 'spans')] = {}
    pairDict[('unique', 'spans')] = {}
    # consider 129,65(F,F)  97,145(F,R) 113,177(R,R)
    # Can these ever be proper paired? What if short indel?
    # need the minus 1 to be exact.
    if str1 == "+":
        pos1left = int(pos1) - 500
        pos1right = int(pos1) + int(repright) - 1
    if str2 == "-":
        pos2left = int(pos2) - 500
        pos2right = int(pos2) + int(repright) - 1
    if str1 == "-":
        pos1left = int(pos1) - int(repright) - 1
        pos1right = int(pos1) + 500
    if str2 == "+":
        pos2left = int(pos2) - int(repright) - 1
        pos2right = int(pos2) + 500
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
                        pairDict[('all', 'spans')][read.query_name] = read.flag
                    if not read.is_duplicate and read.is_paired:
                        pairDict[('unique', 'spans')][read.query_name] = read.flag
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
        pos1left = int(pos1) - 1000
        pos1right = int(pos1) + int(repright) - 1
    if str2 == "-":
        pos2left = int(pos2) - 1000
        pos2right = int(pos2) + int(repright) - 1
    if str1 == "-":
        pos1left = int(pos1) - int(repright) - 1
        pos1right = int(pos1) + 1000
    if str2 == "+":
        pos2left = int(pos2) - int(repright) - 1
        pos2right = int(pos2) + 1000
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


def get_reads_from_bam(bam_file, jxn):
    '''
    Fetch reads from each direction of jxn and get an accounting of each if unique or duplicate.
    '''
    logger.info("Extracting supporting read information for " + jxn)
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
    logger.info("Extracting paired spanning reads for " + jxn)
    # spans
    bam = bamObject.fetch(chrom1, int(pos1) - int(dist1), int(pos1) + dist1)
    spanD = find_discspan(jxn, bam)
    results['spans'] = spanD
    # left junctions
    logger.info("Extracting junction reads originating from the first breakpoint for " + jxn)
    bam = bamObject.fetch(chrom1, int(pos1) - int(dist1), int(pos1) + dist1)
    jxnD_for = find_junctions(bam, chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright)
    results['jxnleft'] = jxnD_for
    # right junctions
    logger.info("Extracting junction reads originating from the second breakpoint for " + jxn)
    flipstr = string.maketrans("-+", "+-")
    nstr1 = str1.translate(flipstr)
    nstr2 = str2.translate(flipstr)
    bam = bamObject.fetch(chrom2, int(pos2) - int(dist2), int(pos2) + dist2)
    jxnD_rev = find_junctions(bam, chrom2, pos2, nstr2, chrom1, pos1, nstr1, repleft, repright)
    results['jxnright'] = jxnD_rev
    bamObject.close()
    logger.info("Finished extracting supporting read information for " + jxn)
    return results


def make_new_dir(newdir):
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass


def subset_bam_by_reads(bam, out_bam, read_ids, cfg):
    logger.info("Subset bam with supporting reads to " + out_bam)
    biobambam2 = (cfg['biobambam2_bin'] + '/bamfilternames')
    names = "names=" + read_ids
    index = "index=1"
    indexfilename = "indexfilename=" + out_bam + ".bai"
    tmpfile = "tmpfile=" + out_bam + ".tmp"
    stdinfile = bam
    bamf = open(out_bam, "wb")
    try:
        with open(stdinfile, "rb") as f:
            retcode = sp.call([biobambam2, names, index, indexfilename,
                               tmpfile], stdin=f, stdout=bamf, stderr=None)
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
    logger.info("Converting bam to fastq for " + in_bam)
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
    logger.info("*velveth Command: " + " ".join(vh_args))
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
    logger.info("*velvetg Command: " + " ".join(vg_args))
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


def do_spades(cfg, assemdir, pfastq, jxnfastq, errlog):
    logger.debug("*Running SPADES")
    spades_cmd = ["python", cfg['spades_exec'], "--12", pfastq, "-s", jxnfastq,
                  "-o", assemdir, "--careful", "--phred-offset", 33,
                  "-t", cfg['spades_threads'], "-m", cfg['spades_mem_gb']]
    spades_args = map(str, spades_cmd)
    logger.info("*SPADES Command: " + " ".join(spades_args))
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


def run_assembly_fxn(jxn, in_bam, cfg, args, *opts):
    ''' Run command to run assembly for a single breakpoint '''
    logger.info("Running assembly fxn in parallel for " + jxn)
    results = {}
    clean_jxn = str(jxn).replace(":", "_")
    clean_jxn = str(clean_jxn).replace("+", "pos")
    clean_jxn = str(clean_jxn).replace("-", "neg")
    jxn_dir = "support" + "/" + clean_jxn + "/"
    make_new_dir(jxn_dir)
    # get supporting reads
    extracted = get_reads_from_bam(in_bam, jxn)
    logger.info("Found and extracted reads successfully")
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
    support_bam = jxn_dir + "/supporting.bam"
    logger.info("*Subsetting reads from bam")
    subset_bam_by_reads(in_bam, support_bam, read_ids, cfg)
    pysam.index(support_bam)
    # convert to fastq and run velvet
    pairfq = jxn_dir + 'paired.fastq'
    junctionfq = jxn_dir + 'junctions.fastq'
    bam2fastq(jxn_dir, support_bam, junctionfq, pairfq)
    errlog = open(jxn_dir + "assembly_log.txt", "w")
    # velvet_all = do_velvet(jxn_dir + "assem_pair", junctionfq, 17, errlog, pairfq)
    velvet_jxn = do_velvet(jxn_dir + "assem_jxn", junctionfq, 17, errlog)
    errlog.close()
    if velvet_jxn:
        vname, vseq = velvet_jxn[0]
        results['velvet'] = vseq
    if args.spades:
        splog = open(jxn_dir + "spades_log.txt", "w")
        spades_seq = do_spades(cfg, jxn_dir + "spades", pairfq, junctionfq, splog)
        splog.close()
        if spades_seq:
            sname, sseq = spades_seq[0]
            results['spades'] = sseq
    results['name'] = jxn
    return results


def init_worker():
    '''this is to sidestep multiproc bug with keyboard interrupts'''
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def run_assembly_parallel(df, in_bam, cfg, args):
    '''Run assembly in parallel'''
    logger.info("Starting Assembly Process")
    make_new_dir("support")
    results = []
    pool = mp.Pool(int(args.workers), init_worker)
    try:
        for jxn in df['name']:
            seq = pool.apply_async(run_assembly_fxn, args=[jxn, in_bam, cfg, args])
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
    logger.info("Finished Assembly")
    return list_of_dicts
