#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import re
import logging
import subprocess as sp
import pysam  # requires 0.9.0 or newer

__author__ = "Jeff Jasper"
__email__ = "jasper1918@gmail.com"

logger = logging.getLogger("STAR-SEQR")


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
            STAR_args = ['STAR', '--readFilesIn', fq1, fq2, '--readFilesCommand', 'zcat',
                         '--runThreadN', args.threads, '--genomeDir', cfg['star_index_dna'],
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
                sens_params = ['--chimSegmentMin', 7, '--chimJunctionOverhangMin', 7,
                               '--chimScoreMin', 1, '--chimScoreDropMax', 20,
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
            STAR_args = ['STAR', '--readFilesIn', fq1, fq2, '--readFilesCommand', 'zcat',
                         '--runThreadN', args.threads, '--genomeDir', cfg['star_index_rna'],
                         '--outFileNamePrefix ', args.prefix + ".", '--outSAMtype', 'None',
                         '--alignIntronMax', 200000, '--alignMatesGapMax', 200000,
                         '--chimOutType', 'SeparateSAMold', '--chimScoreJunctionNonGTAG', -1,
                         '--alignSJDBoverhangMin', 10]  # , '--outSJfilterCountTotalMin', '5 -1 5 5']
            # choose sensitivity mode
            if (args.mode == 0):
                sens_params = ['--chimSegmentMin', 12, '--chimJunctionOverhangMin', 12,
                               '--chimScoreMin', 1, '--chimScoreDropMax', 20,
                               '--chimScoreSeparation', 10, '--chimSegmentReadGapMax', 3,
                               '--chimFilter', 'None', '--twopassMode', "Basic"]
            elif (args.mode == 1):
                sens_params = ['--chimSegmentMin', 5, '--chimJunctionOverhangMin', 5,
                               '--chimScoreMin', 0, '--chimScoreDropMax', 40,
                               '--chimScoreSeparation', 10, '--chimSegmentReadGapMax', 3,
                               '--chimFilter', 'None', '--twopassMode', "Basic"]
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
    mark_samdups_cmd = ['bammarkduplicates2', "=".join(["I", in_bam]), "=".join(["O", out_bam]), "=".join(["markthreads",
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
