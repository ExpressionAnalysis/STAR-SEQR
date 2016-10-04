#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import time
import re
import errno
import string
from argparse import ArgumentParser
import pandas as pd
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import gzip
import io
import logging
import numpy as np
import multiprocessing as mp
import signal
import starseqr_utils


def parse_args():
    usage = " "
    parser = ArgumentParser(description="STAR-SEQR Parameters:", epilog=usage)
    # create STAR alignment
    group1 = parser.add_argument_group('Do Alignment', '')
    group1.add_argument('-1', '--fastq1', type=str, required=False,
                        help='fastq 1')
    group1.add_argument('-2', '--fastq2', type=str, required=False,
                        help='fastq 2')
    group1.add_argument('-i', '--star_index', type=str, required=False,
                        help='path to STAR index folder')
    group1.add_argument('-m', '--mode', type=int, required=False,
                        default=0,
                        choices=[0, 1],
                        help='STAR alignment sensitivity Mode. 0=Default, 1=More-Sensitive')
    # existing STAR alignment
    group2 = parser.add_argument_group('Use Existing Alignment', '')
    group2.add_argument('-sj', '--star_jxns', type=str, required=False,
                        help='chimeric junctions file produce by STAR')
    group2.add_argument('-ss', '--star_sam', type=str, required=False,
                        help='chimeric sam file produced by STAR')
    # shared args
    parser.add_argument('-p', '--prefix', type=str, required=True,
                        help='prefix to name files')
    parser.add_argument('-d', '--dist', type=int, required=False,
                        default=100000,
                        help='minimum distance to call junctions')
    parser.add_argument('-j', '--jxn_reads', type=int, required=False,
                        default=2,
                        help='minimum number of junction reads to keep junctions.')
    parser.add_argument('-s', '--span_reads', type=int, required=False,
                        default=2,
                        help='minimum number of spanning discordant read pairs to call junctions.')
    parser.add_argument('-n', '--nucleic_type', type=str, required=False,
                        default="RNA",
                        help='nucleic acid type',
                        choices=["RNA", "DNA"])
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=12,
                        help='Number of threads to use for STAR')
    parser.add_argument('--ann_source', '--ann_source', type=str, required=False,
                        default="gencode",
                        help='annotation source',
                        choices=["refgene", "ensgene", "gencode"])
    parser.add_argument('-b', '--bed_file', type=str, required=False,
                        help='Bed file to subset analysis')
    parser.add_argument('--spades', '--spades', action='store_true',
                        help='do spades assembly instead of velvet')
    parser.add_argument('--keep_dups', action='store_true',
                        help='keep read duplicates')
    parser.add_argument('--keep_gene_dups', action='store_true',
                        help='allow RNA internal gene duplications to be considered')
    parser.add_argument('--keep_mito', action='store_true',
                        help='allow RNA fusions to contain at least one breakpoint from Mitochondria')
    parser.add_argument('--keep_novel', action='store_true',
                        help='allow gene fusions to pass that have no gene annotation')
    parser.add_argument('--keep_unscaffolded', action='store_true',
                        help='allow gene fusions to pass that are found on unscaffolded contigs')
    parser.add_argument('-v', '--verbose', action="count",
                        help="verbose level... repeat up to three times.")
    args = parser.parse_args()

    # check that the correct args have been specified
    align = [args.fastq1, args.fastq2, args.star_index]
    call = [args.star_jxns, args.star_sam]
    if any(align) and any(call):
        print("Error: Please choose either fastqs or STAR existing files as input!")
        sys.exit(1)
    if any(align) and None in align:
        print("Error: Fastq1, Fastq2, and the STAR index must be specified if doing alignment")
        sys.exit(1)
    if any(call) and None in call:
        print("Error: The STAR .junctions and .sam file must be specified if using existing alignment")
        sys.exit(1)
    return args


def check_file_exists(path):
    if (os.stat(os.path.realpath(path)).st_size == 0):
        logger.error("Exiting. Cannot find file: " + os.path.realpath(path))
        sys.exit(1)

def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)


def find_resource(filename):
    packagedir = starseqr_utils.__path__[0]
    dirname = os.path.join(packagedir, 'resources')
    fullname = os.path.abspath(os.path.join(dirname, filename))
    return fullname


def set_log_level_from_verbose(ch, args):
    logger = logging.getLogger("STAR-SEQR")
    if not args.verbose:
        ch.setLevel('ERROR')
    elif args.verbose == 1:
        ch.setLevel('WARNING')
    elif args.verbose == 2:
        ch.setLevel('INFO')
    elif args.verbose >= 3:
        ch.setLevel('DEBUG')
    else:
        logger.critical("UNEXPLAINED VERBOSE VALUE")


def import_jxns_pandas(jxnFile, args):
    logger = logging.getLogger("STAR-SEQR")
    logger.info('Importing junctions')
    df = pd.read_csv(jxnFile, sep="\t", header=None, usecols=range(0, 14), low_memory=False,  engine='c')
    df.columns = ['chrom1', 'pos1', 'str1', 'chrom2', 'pos2', 'str2', 'jxntype', 'jxnleft', 'jxnright', 'readid', 'base1', 'cigar1', 'base2', 'cigar2']
    df['readid'] = df['readid'].astype(str)
    df['pos1'] = df['pos1'].astype(int)
    df['pos2'] = df['pos2'].astype(int)
    df['identity'] = df['base1'].astype(str) + ':' + df['cigar1'].astype(str) + ':' + df['base2'].astype(str) + ':' + df['cigar2'].astype(str)
    df.drop(['base1', 'cigar1', 'base2', 'cigar2'], axis=1, inplace=True)
    if args.keep_dups:
        logger.info("Allowing duplicate reads")
        return df
    else:
        logger.info("Removing duplicate reads")
        return df.drop_duplicates(subset=['identity'], keep='first')


def pandas_parallel(df, func, nthreads):
    logger = logging.getLogger("STAR-SEQR")
    start = time.time()
    def init_worker():
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        df_split = np.array_split(df, min(nthreads, len(df.index)))
        pool = mp.Pool(nthreads, init_worker)
        these_res = pool.map(func, df_split)
        df = pd.concat(these_res)
        pool.close()
        pool.join()
    except KeyboardInterrupt as e:
        logger.error("Error: Keyboard interrupt")
        pool.terminate()
        raise e
    except Exception as e:
        logger.error("Exception: " + str(e))
        pool.terminate()
        raise e
    logger.info("Time to run pandas_parallel on " + str(func.__name__) + " took %g seconds" %  (time.time() - start))
    return df


def apply_choose_order(df):
    df['order'] = df.apply(lambda x: choose_order(x['chrom1'], x['pos1'], x['chrom2'], x['pos2']), axis=1)
    return df


def choose_order(chrL1, posL1, chrL2, posL2):
    ''' Choose one reprentation of jxn to merge on.
    This is just for DNA breakpoints where bidirectional breakpoint detection occurs'''
    mychrL1 = str(chrL1).replace("chr", "")
    mychrL2 = str(chrL2).replace("chr", "")
    chrList = [mychrL1, mychrL2]
    chrList.sort()
    if ((mychrL1) == (mychrL2)):
        if (int(posL1) < int(posL2)):
            return 1
        else:
            return 2
    elif (int(chrList.index(mychrL1)) < int(chrList.index(mychrL2))):
        return 1
    elif (int(chrList.index(mychrL1)) > int(chrList.index(mychrL2))):
        return 2


def apply_normalize_jxns(df):
    df['name'] = df.apply(lambda x: normalize_jxns(x['chrom1'], x['chrom2'], x['pos1'], x['pos2'],
                                                           x['str1'], x['str2'], x['jxnleft'], x['jxnright'],
                                                           x['order']), axis=1)
    return df


def normalize_jxns(chrom1, chrom2, pos1, pos2, strand1, strand2, repleft, repright, order):
    '''Choose one representation for DNA breakpoints'''
    flipstr = string.maketrans("-+", "+-")
    if order == 2:
        if strand1 == "-":
            new_pos1 = str(chrom1) + ":" + str(pos1 - int(repright)) + ":" + strand1.translate(flipstr)
        else:
            new_pos1 = str(chrom1) + ":" + str(pos1 + int(repright)) + ":" + strand1.translate(flipstr)
        if strand2 == "-":
            new_pos2 = str(chrom2) + ":" + str(pos2 - int(repright)) + ":" + strand2.translate(flipstr)
        else:
            new_pos2 = str(chrom2) + ":" + str(pos2 + int(repright)) + ":" + strand2.translate(flipstr)
        newid = new_pos2 + ":" + new_pos1 + ":" + str(repleft) + ":" + str(repright)
    elif order == 1:
        new_pos1 = str(chrom1) + ":" + str(pos1) + ":" + strand1
        new_pos2 = str(chrom2) + ":" + str(pos2) + ":" + strand2
        newid = new_pos1 + ":" + new_pos2 + ":" + str(repleft) + ":" + str(repright)
    return newid


def parallel_count_jxns(df, args):
    # logger = logging.getLogger("STAR-SEQR")
    # start = time.time()
    # def init_worker():
    #     signal.signal(signal.SIGINT, signal.SIG_IGN)
    # def parallel_groupby(dfGrouped, func, nthreads):
    #     try:
    #         pool = mp.Pool(nthreads, init_worker)
    #         ret_list = pool.map(func, [group for name, group in dfGrouped])
    #         pool.close()
    #         pool.join()
    #         logger.info("Time to aggregate junctions took %g seconds" %  (time.time() - start))
    #         tmp = pd.concat(ret_list)
    #     except KeyboardInterrupt as e:
    #         logger.error("Error: Keyboard interrupt")
    #         pool.terminate()
    #         raise e
    #     except Exception as e:
    #         logger.error("Exception: " + str(e))
    #         pool.terminate()
    #         raise e

    grouped_df = df.groupby(['name', 'order'], as_index=True)
    # new_df = parallel_groupby(grouped_df, apply_count_vals, args.threads)
    new_df = grouped_df.agg({'readid': {'reads': lambda col: ','.join(col), 'counts' :'count'}}).reset_index().pivot(index='name', columns='order').reset_index()
    if args.nucleic_type == "DNA":
        new_df.columns = ['name', 'jxnreadsleft', 'jxnreadsright', 'jxnleft', 'jxnright']
    elif args.nucleic_type == "RNA":
        new_df.columns = ['name', 'jxn_reads', 'jxn_counts']
    return new_df


def apply_count_vals(df):
    # print(df.head(10))
    new_df = df.groupby(['name', 'order'], as_index=True).agg({'readid': {'reads': lambda col: ','.join(col), 'counts' :'count'}}).reset_index().pivot(index='name', columns='order').reset_index()
    return new_df


def bed_to_tree(bed):
    with open(bed, 'r') as f:
        gtree = GenomeIntervalTree.from_bed(fileobj=f)
    return gtree


def subset_bed_func(jxn, targets_tree):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    intersect1 = len(targets_tree[str(chrom1)].search(int(pos1)))
    intersect2 = len(targets_tree[str(chrom2)].search(int(pos2)))
    if intersect1 or intersect2 >= 1:
        return 1
    else:
        return 0


def get_distance(jxn):
    ''' reference distance'''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if chrom1 == chrom2:
        dist_between = abs(int(pos1) - int(pos2))
    else:
        dist_between = np.nan  # max value to use between chromosomes.
    return dist_between


def apply_pairs_func(df):
    df['spans'], df['spanreads'] = zip(*df.apply(lambda x: get_pairs_func(x['name']), axis=1))
    return df


def apply_primers_func(df):
    df['primers'] = df.apply(lambda x: starseqr_utils.run_primer3.runp3(x['name'], x['assembly']), axis=1).apply(lambda x: ",".join(x))
    return df['primers']


def get_pairs_func(jxn):
    '''
    Get paired end read data that supports each jxn from the junction file.
    These have a -1 for jxntype.
    Need to grab both sides of jxn. Each is unique in the .junction file.
    Querying the file can be tedious with the flipping of the junctions.
    '''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1 = int(pos1)
    pos2 = int(pos2)
    maxhom = max(int(repright), int(repleft))
    if str1 == "+":
        pos1left = pos1 - 100000
        pos1right = pos1 + maxhom + 1
    elif str1 == "-":
        pos1left = pos1 - maxhom - 1
        pos1right = pos1 + 100000
    if str2 == "+":
        pos2left = pos2 - maxhom - 1
        pos2right = pos2 + 100000
    elif str2 == "-":
        pos2left = pos2 - 100000
        pos2right = pos2 + maxhom + 1
    forward = rawdf[(rawdf['jxntype'] == -1) &
                    (rawdf['chrom1'] == chrom1) & (rawdf['chrom2'] == chrom2) &
                    (rawdf['pos1'] >= pos1left) & (rawdf['pos1'] <= pos1right) &
                    (rawdf['pos2'] >= pos2left) & (rawdf['pos2'] <= pos2right)]
    reverse = rawdf[(rawdf['jxntype'] == -1) &
                    (rawdf['chrom1'] == chrom2) & (rawdf['chrom2'] == chrom1) &
                    (rawdf['pos1'] >= pos2left) & (rawdf['pos1'] <= pos2right) &
                    (rawdf['pos2'] >= pos1left) & (rawdf['pos2'] <= pos1right)]
    npairs = len(forward['readid'].index) + len(reverse['readid'].index)
    reads = ','.join(forward['readid'].tolist()) + ',' +  ','.join(reverse['readid'].tolist())
    return (npairs, reads)


def apply_jxn_strand(df):
    _, _, df['test_strand'] = zip(*df.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
    return df


def apply_flip_func(df):
    df['name'], df['flip'] = zip(*df.apply(lambda x: flip_jxn(x['name'], x['test_strand']), axis=1))
    return df


def flip_jxn(jxn, gs1):
    '''Flip jxn orientation for RNA Fusion that is inverse according to strand info'''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if str1 != gs1[0]:
        flip = 1
        flipstr = string.maketrans("-+", "+-")
        if str1 == "-":
            new_pos1 = str(chrom1) + ":" + str(int(pos1)) + ":" + str1.translate(flipstr)
        else:
            new_pos1 = str(chrom1) + ":" + str(int(pos1)) + ":" + str1.translate(flipstr)
        if str2 == "-":
            new_pos2 = str(chrom2) + ":" + str(int(pos2)) + ":" + str2.translate(flipstr)
        else:
            new_pos2 = str(chrom2) + ":" + str(int(pos2)) + ":" + str2.translate(flipstr)
        newid = new_pos2 + ":" + new_pos1 + ":" + str(repright) + ":" + str(repleft)
    else:
        newid = jxn
        flip = 0
    return (newid, flip)


def get_svtype_func(jxn):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if chrom1 != chrom2:
        svtype = "TRANSLOCATION"
    else:
        # STAR notation is same as other tools after strand2 is flipped.
        if str(str1) == "+" and str(str2) == "-":
            svtype = "INSERTION"
        elif str(str1) == "-" and str(str2) == "+":
            svtype = "INVERSION"
        elif str(str1) == "+" and str(str2) == "+":
            svtype = "DELETION"
        elif str(str1) == "-" and str(str2) == "-":
            svtype = "DUPLICATION"
    return svtype


def get_sv_locations(jxn):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1 = int(pos1) - int(repright)
    pos2 = int(pos2) # + repleft?
    if str(str1) == "+" and str(str2) == "+":
        pos1 -= 1
        pos2 += 1
    elif str(str1) == "-" and str(str2) == "-":
        pos1 += 1
        pos2 -= 1
    elif str(str1) == "+" and str(str2) == "-":
        pos1 -= 1
        pos2 -= 1
    elif str(str1) == "-" and str(str2) == "+":
        pos1 += 1
        pos2 += 1
    brk1 = str(chrom1) + ":" + str(pos1)
    brk2 = str(chrom2) + ":" + str(pos2)
    return (brk1, brk2)


def main():
    start = time.time()
    args = parse_args()

    # file log
    cust_format = '%(asctime)s - %(module)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=cust_format,
                        datefmt="%Y-%m-%d %H:%M:%S", filename=args.prefix + '_STAR-SEQR.log')
    # console log
    logger = logging.getLogger("STAR-SEQR")
    console = logging.StreamHandler()
    set_log_level_from_verbose(console, args)
    formatter = logging.Formatter('%(asctime)-15s - %(levelname)-4s - %(message)s', "%Y-%m-%d %H:%M")
    console.setFormatter(formatter)
    logger.addHandler(console)

    # start analysis
    logger.info("***************STAR-SEQR******************")
    logger.info("CMD = " + str(' '.join(sys.argv)))
    logger.info('Starting to work on sample: ' + args.prefix)

    # check files exist and get abs paths
    if args.bed_file:
        check_file_exists(bed_file)
        bed_path = os.path.realpath(args.bed_file)
    if args.fastq1:
        fq1_path = os.path.realpath(args.fastq1)
        fq2_path = os.path.realpath(args.fastq2)
        check_file_exists(fq1_path)
        check_file_exists(fq2_path)
    if args.star_jxns:
        starjxns_path = os.path.realpath(args.star_jxns)
        starsam_path = os.path.realpath(args.star_sam)
        check_file_exists(starjxns_path)
        check_file_exists(starsam_path)

    # make sample folder
    if not os.path.exists(args.prefix + "_STAR-SEQR"):
        os.makedirs(args.prefix + "_STAR-SEQR")
    os.chdir(args.prefix + "_STAR-SEQR")

    # Do alignment if fastqs
    if args.fastq1:
        starseqr_utils.star_funcs.run_star(fq1_path, fq2_path, args)
    # symlink existing files into folder if provided
    elif args.star_jxns:
        force_symlink(starjxns_path, args.prefix + ".Chimeric.out.junction")
        force_symlink(starsam_path, args.prefix + ".Chimeric.out.sam")

    # import all jxns
    global rawdf
    rawdf = import_jxns_pandas(args.prefix + ".Chimeric.out.junction", args)
    jxns = rawdf[rawdf['jxntype'] >= 0].reset_index() # junctions can be either 0, 1, 2

    if len(jxns.index) == 0:
        logger.info("No junctions found in the input file")
        starseqr_utils.sv2bedpe.process(jxns, args)
        sys.exit(0)

    # Prepare Annotation
    global gtree
    if args.ann_source == "refgene":
        refgene = find_resource("refGene.txt.gz")
        kg = io.BufferedReader(gzip.open(refgene))
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.REF_GENE)
    elif args.ann_source == "ensgene":
        ensgene = find_resource("ensGene.txt.gz")
        kg = io.BufferedReader(gzip.open(ensgene))
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)
    elif args.ann_source == "gencode":
        gencode = find_resource("wgEncodeGencodeBasicV24lift37.txt.gz")
        kg = io.BufferedReader(gzip.open(gencode))
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)

    if args.nucleic_type == "RNA":
        # start output files
        stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')
        breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
        breakpoint_cols = ["ann", "breakpoint_left", "breakpoint_right", "left_symbol", "right_symbol", "left_annot", "right_annot", "dist", "spans_disc_all", "jxn_first", "jxn_second", "assembly", "primers", "name"]
        breakpoint_header = ["NAME", "BRKPT_LEFT", "BRKPT_RIGHT", "LEFT_SYMBOL", "RIGHT_SYMBOL", "LEFT_ANNOT", "RIGHT_ANNOT", "DISTANCE", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT", "ASSEMBLY", "PRIMERS", "ID"]
        print(*breakpoint_header, sep='\t', file=breakpoints_fh)

        # stats dict
        stats_res = {'Total_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}

         # Order, Normalize and Aggregate
        logger.info("Ordering junctions")
        jxns['order'] = 1
        logger.info('Normalizing junctions')
        jxns = pandas_parallel(jxns, apply_normalize_jxns, args.threads)
        logger.info("Getting gene strand and flipping info as necessary")
        jxns = pandas_parallel(jxns, apply_jxn_strand, args.threads)
        jxns = pandas_parallel(jxns, apply_flip_func, args.threads)
        logger.info("Aggregating junctions")
        jxn_summary = parallel_count_jxns(jxns, args)

        # write stats
        stats_res['Total_Breakpoints'] = len(jxn_summary.index)
        logger.info('Total Breakpoints:' + str(stats_res['Total_Breakpoints']))

        # Get discordant pairs
        logger.info('Getting pair info')
        jxn_filt = pandas_parallel(jxn_summary, apply_pairs_func, args.threads)

        # hard filter on minimal junction and span reads
        logger.info('Filtering junctions')
        jxn_filt = jxn_filt[(jxn_filt["spans"] *  3 + jxn_filt["jxn_counts"] * 3)  >= 6]

        if len(jxn_filt.index) >= 1:
            # Get dist
            jxn_filt['dist'] = jxn_filt.apply(lambda x: get_distance(x['name']), axis=1)

            # Get Annotation info for each junction
            jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'] = zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
            jxn_filt['right_symbol'], jxn_filt['right_annot'], jxn_filt['right_strand']= zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_anno(x['name'], gtree, 2), axis=1))
            # get all genes associated to look for overlap for each read later..
            jxn_filt['txinfo'] = jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxn_info_func(x['name'], gtree), axis=1)

            # subset to ROI using bed file if it exists
            if args.bed_file:
                logger.info('Subsetting junctions using the supplied bed file')
                targets_tree = bed_to_tree(bed_path)
                jxn_filt['subset'] = jxn_filt.apply(lambda x: subset_bed_func(x['name'], targets_tree), axis=1)
                jxn_filt = jxn_summary[jxn_filt['subset'] >= 1]

            # remove internal gene dups unless otherwise requested
            if not args.keep_gene_dups:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[jxn_filt['left_symbol'] != jxn_filt['right_symbol']]  # todo: this should rather do a set overlap of alltranscripts found.
                logger.info("Number of candidates removed due to internal gene duplication filter: " + str(before_remove - len(jxn_filt.index)))

            # remove novel genes unless otherwise requested
            if not args.keep_novel:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[(jxn_filt['left_symbol'] != "NA") & (jxn_filt['right_symbol'] != "NA")]  # todo: this should rather do a set overlap of alltranscripts found.
                logger.info("Number of candidates removed due to novel gene filter: " + str(before_remove - len(jxn_filt.index)))

            # remove mitochondria unless otherwise requested
            if not args.keep_mito:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[~jxn_filt['name'].str.contains("chrM|M")]
                logger.info("Number of candidates removed due to Mitochondria filter: " + str(before_remove - len(jxn_filt.index)))

            # remove non-canonical contigs unless otherwise requested
            if not args.keep_unscaffolded:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[~jxn_filt['name'].str.contains("Un|random|hap")]
                logger.info("Number of candidates removed due to non-canonical contigs filter: " + str(before_remove - len(jxn_filt.index)))

            #combine all supporting reads together.
            jxn_filt['supporting_reads'] = jxn_filt['spanreads']  + ',' + jxn_filt['jxn_reads']

            stats_res['Candidate_Breakpoints'] = len(jxn_filt.index)
            logger.info('Candidate Breakpoints:' + str(stats_res['Candidate_Breakpoints']))

        # Process candidates
        if len(jxn_filt.index) >= 1:
            # identify duplicate reads and mark the chimeric fragment
            starseqr_utils.star_funcs.convert(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.bam", args)
            # Gather unique read support
            assemdict = starseqr_utils.support_funcs_rna.run_support_parallel(jxn_filt, args.prefix + ".Chimeric.out.bam", args)
            logger.info("Finished aggregating support.")
            assemdf = pd.DataFrame.from_records(assemdict, index='name')
            # Merge read support with previous stats
            finaldf = pd.merge(jxn_filt, assemdf, how='inner', left_on="name", right_on="name", left_index=False,
                               right_index=True, sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
            # collapse read info for brevity but keep here in case useful later on
            finaldf['jxn_first'] = finaldf["jxnleft_for_first"] + finaldf["jxnleft_rev_first"] + \
                                   finaldf["jxnright_for_first"] + finaldf["jxnright_rev_first"]
            finaldf['jxn_second'] = finaldf["jxnleft_for_second"] + finaldf["jxnleft_rev_second"] + \
                                    finaldf["jxnright_for_second"] + finaldf["jxnright_rev_second"]
            finaldf['spans_disc'] = finaldf['spans_disc_all']
            # todo: confirm breakpoint with bwa or bowtie or age?
            # todo: get novel or existing fusion info from Fusco.
            # todo: probabilistic module to assign quality score

            # Generate Primers using assembled contigs
            logger.info("Generating primers from assembled contigs")
            finaldf['primers'] = pandas_parallel(finaldf, apply_primers_func, args.threads)

            # Get Annotation
            finaldf['ann'] = finaldf['left_symbol'] + "--" + finaldf['right_symbol']

            # Get breakpoint locations
            finaldf['breakpoint_left'], finaldf['breakpoint_right']  = zip(*finaldf.apply(lambda x: get_sv_locations(x['name']), axis=1))

            # all candidates
            finaldf.to_csv(path_or_buf="STAR-SEQR_candidate_info.txt", header=True, sep="\t", mode='w', index=False)

            # Hard filter on read counts after accounting for transcript info. Change this once a probabilistic module is ready.
            finaldf = finaldf[(finaldf["spans_disc_all"] *  4 + finaldf["jxn_first"] * 2 + finaldf["jxn_second"] * 2)  >= 6]

            # Write output
            finaldf.sort_values(['jxn_first', "spans_disc_all"], ascending=[False, False], inplace=True)
            finaldf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                           columns=breakpoint_cols, mode='w', index=False)
            breakpoints_fh.close()

            # Make bedpe and VCF
            starseqr_utils.sv2bedpe.process(finaldf, args)

            # Log Stats
            stats_res['Passing_Breakpoints'] = len(finaldf.index)
            logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))

        if len(jxn_filt.index) == 0:
            logger.info("No junctions found after filtering")
            starseqr_utils.sv2bedpe.process(jxn_filt, args)
            sys.exit(0)

    elif args.nucleic_type == "DNA":
        # start output files
        stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')
        breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
        breakpoint_cols = ["ann", "svtype", "breakpoint_left", "breakpoint_right", "dist", "spans_disc", "jxn_first", "jxn_second", "name"]
        breakpoint_header = ["NAME", "SVTYPE", "BRKPT_LEFT", "BRKPT_RIGHT", "DISTANCE", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT", "ID"]
        print(*breakpoint_header, sep='\t', file=breakpoints_fh)

        # stats dict
        stats_res = {'Total_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}

        # Order, Normalize and Aggregate
        logger.info("Ordering junctions")
        jxns = pandas_parallel(jxns, apply_choose_order, args.threads)
        logger.info('Normalizing junctions')
        jxns = pandas_parallel(jxns, apply_normalize_jxns, args.threads)
        logger.info("Aggregating junctions")
        jxn_summary = parallel_count_jxns(jxns, args)

        # subset to bed if specified
        if args.bed_file:
            start_bed = time.time()
            logger.info('Subsetting junctions using the supplied bed file')
            targets_tree = bed_to_tree(bed_path)
            jxn_summary['subset'] = jxn_summary.apply(lambda x: subset_bed_func(x['name'], targets_tree), axis=1)
            jxn_summary = jxn_summary[jxn_summary['subset'] >= 1]
            logger.info("Time to subset junction from bed took %g seconds" %  (time.time() - start_bed))

        # print stats
        stats_res['Total_Breakpoints'] = len(jxn_summary.index)
        logger.info('Total Breakpoints:' + str(stats_res['Total_Breakpoints']))

        # Filter on distance
        logger.info('Filtering junctions based on distance')
        jxn_summary['dist'] = jxn_summary.apply(lambda x: get_distance(x['name']), axis=1)
        jxn_summary = jxn_summary[(jxn_summary['dist'] >= args.dist) |
                                  (pd.isnull(jxn_summary['dist']))]
        logger.info('Junctions passing distance filter:' + str(len(jxn_summary.index)))

        # Filter on junctions
        logger.info('Filtering junctions based on number of junction reads')
        jxn_filt = jxn_summary[(jxn_summary["jxnleft"] >= args.jxn_reads) &
                            (jxn_summary["jxnright"] >= args.jxn_reads)].sort_values("jxnleft", ascending=False)
        logger.info('Junctions passing junction read filter:' + str(len(jxn_filt.index)))


        if len(jxn_filt.index) >= 1:
             # Annotate genes
            jxn_filt['genesleft'], jxn_filt['genesright'], jxn_filt['common'] = zip(*jxn_filt.apply(lambda x: starseqr_utils.annotate_sv.get_jxnside_genes(x['name'], gtree), axis=1))

            # Get gene info and remove internal gene dups
            if not args.keep_gene_dups:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[jxn_filt['common'] == 0]
                logger.info("Number of candidates removed due to internal gene duplication filter: " + str(before_remove - len(jxn_filt.index)))


        if len(jxn_filt.index) >= 1:
            # Get paired discordant spanning reads supporting junctions and filter
            logger.info('Getting pair info')
            jxn_filt = pandas_parallel(jxn_filt, apply_pairs_func, args.threads)
            logger.info('Filtering junctions based on pairs')
            jxn_filt = jxn_filt[(jxn_filt['spans'] >= args.span_reads)]

            #combine all supporting reads together.
            jxn_filt['supporting_reads'] = jxn_filt['spanreads']  + ',' + jxn_filt['jxnreadsleft']   + ',' + jxn_filt['jxnreadsright']

            # Write candidates to file
            # jxn_filt.to_csv(path_or_buf="STAR-SEQR_candidates.txt", header=True, sep="\t", mode='w',
            #                  index=False, columns=["name", "jxnleft", "jxnright", "spans", "dist", "supporting_reads"])

            # Log Stats
            stats_res['Candidate_Breakpoints'] = len(jxn_filt.index)
            logger.info('Candidate Breakpoints:' + str(stats_res['Candidate_Breakpoints']))

        if len(jxn_filt.index) >= 1:
            # skip the big functions for now
            finaldf = jxn_filt
            finaldf['jxn_first'] = finaldf['jxnleft']
            finaldf['jxn_second'] = finaldf['jxnright']
            finaldf['spans_disc'] = finaldf['spans']

            # mark duplicates
            starseqr_utils.star.convert(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.bam", args)

            # # Get read support from BAM
            # assemdict = assem_dna.run_support_parallel(jxn_filt, args.prefix + ".Chimeric.out.bam", args)
            # assemdf = pd.DataFrame.from_records(assemdict, index='name')

            # # Merge Stats for use
            # finaldf = pd.merge(jxn_filt, assemdf, how='inner', left_on="name", right_on="name", left_index=False,
            #                    right_index=True, sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
            # finaldf['jxn_first'] = finaldf["jxnleft_for_first"] + finaldf["jxnleft_rev_first"] + \
            #     finaldf["jxnright_for_first"] + finaldf["jxnright_rev_first"]
            # finaldf['jxn_second'] = finaldf["jxnleft_for_second"] + finaldf["jxnleft_rev_second"] + \
            #     finaldf["jxnright_for_second"] + finaldf["jxnright_rev_second"]
            # finaldf['spans_disc'] = finaldf['spans_disc_all']

            # todo: confirm breakpoint
            # todo: probabilist module to assign quality score

            # Generate Primers
            # finaldf['primers'] = pandas_parallel(finaldf, apply_primers_func, args.threads)

            # Get Annotation
            finaldf['ann'] = starseqr_utils.annotate_sv.get_gene_info(finaldf, gtree)

            # Get Breakpoint type
            finaldf['svtype'] = finaldf.apply(lambda x: get_svtype_func(x['name']), axis=1)

            # Get breakpoint locations
            finaldf['breakpoint_left'], finaldf['breakpoint_right']  = zip(*finaldf.apply(lambda x: get_sv_locations(x['name']), axis=1))

            finaldf.to_csv(path_or_buf="STAR-SEQR_candidate_info.txt", header=True, sep="\t", mode='w', index=False)
            # finaldf = finaldf[(finaldf["spans_disc_unique"] >= args.span_reads) &
            #                   (finaldf["jxn_first_unique"] >= args.jxn_reads) &
            #                   (finaldf["jxn_second_unique"] >= args.jxn_reads)]

            # Write output
            finaldf.sort_values(['jxnleft', "spans"], ascending=[False, False], inplace=True)
            finaldf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                           columns=breakpoint_cols, mode='w', index=False)
            breakpoints_fh.close()

            # Make bedpe and VCF
            starseqr_utils.sv2bedpe.process(finaldf, args)

            # Log Stats
            stats_res['Passing_Breakpoints'] = len(finaldf.index)
            logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))
        else:
            logger.info("No candidate junctions identified.")
            # Make bedpe and VCF, write headers only
            starseqr_utils.sv2bedpe.process(jxn_filt, args)


    # Write stats to file
    for key, value in stats_res.items():
        stats_fh.write(key + "\t" + str(value) + "\n")


    # Finish
    logger.info("Program took  %g seconds" % (time.time() - start))


if __name__ == "__main__":
    sys.exit(main())
