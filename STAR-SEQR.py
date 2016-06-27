#!/usr/bin/env python

import os
import sys
import time
import re
import string
from ConfigParser import SafeConfigParser
from argparse import ArgumentParser
import pandas as pd
from intervaltree_bio import GenomeIntervalTree
import aux_funcs as aux
import annotate_sv as ann
import sv2bedpe
import run_primer3 as primer3
import logging

__author__ = "Jeff Jasper"
__email__ = "jasper1918@gmail.com"


def parse_args():
    usage = " "
    parser = ArgumentParser(
        description="STAR-SEQR Parameters:", epilog=usage)
    parser.add_argument('-1', '--fastq1', type=str, required=True,
                        help='fastq 1',
                        metavar="fastq 1")
    parser.add_argument('-2', '--fastq2', type=str, required=False,
                        help='fastq 2',
                        metavar="fastq 2")
    parser.add_argument('-p', '--prefix', type=str, required=True,
                        help='prefix to name files',
                        metavar="prefix")
    parser.add_argument('-d', '--dist', type=int, required=False, default=500,
                        help='minimum distance to call junctions',
                        metavar="distance_threshold")
    parser.add_argument('-j', '--jxn_reads', type=int, required=False,
                        default=2,
                        help='minimum number of junction reads to keep junctions.',
                        metavar="junction_depth")
    parser.add_argument('-s', '--span_reads', type=int, required=False,
                        default=2,
                        help='minimum number of spanning discordant read pairs to call junctions.',
                        metavar="spanning_depth")
    parser.add_argument('-n', '--nucleic_type', type=str, required=False,
                        default="DNA",
                        help='nucleic acid type (DNA, RNA)',
                        metavar="nucleic acid type")
    parser.add_argument('-m', '--mode', type=int, required=False,
                        default=0,
                        help='Sensitivity Mode. 0=Default, 1=Extra-Sensitive',
                        metavar="STAR Parameters to invoke to make more sensitive.")
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=12,
                        help='Number of threads to use for STAR',
                        metavar="Threads")
    parser.add_argument('--bidir', '--bidir', action='store_true',
                        help='require bidirectional breakpoints for detection')
    parser.add_argument('-b', '--bed_file', type=str, required=False,
                        help='Bed file to subset analysis',
                        metavar="Bed file, 3 col minimal")
    parser.add_argument('-c', '--config_file', type=str, required=True,
                        help='Config file of paths to key programs and files',
                        metavar="config file.ini")
    parser.add_argument('-w', '--workers', type=int, required=False,
                        default=3,
                        help='number of workers',
                        metavar="number of workers")
    parser.add_argument('--spades', '--spades', action='store_true',
                        help='do spades assembly')
    parser.add_argument('-v', '--verbose', action="count",
                        help="verbose level... repeat up to three times.")
    return parser.parse_args()


def parse_config(config_file):
    logger = logging.getLogger("STAR-SEQR")
    if (os.stat(os.path.realpath(config_file)).st_size == 0):
        logger.critical("Exiting. Cannot find file: " + os.path.realpath(config_file))
        sys.exit(1)
    cfgparser = SafeConfigParser()
    cfgparser.read(config_file)
    # print(json.dumps(cfgparser._sections['PATHS'], sort_keys=True, indent=4))
    return cfgparser._sections['PATHS']


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


def import_jxns_pandas(jxnFile):
    logger = logging.getLogger("STAR-SEQR")
    logger.info('Importing junctions')
    df = pd.read_csv(jxnFile, sep="\t", header=None, usecols=range(0, 10))
    df.columns = ['chrom1', 'pos1', 'str1', 'chrom2', 'pos2', 'str2', 'jxntype', 'jxnleft', 'jxnright', 'readid']
    return df


def choose_order(chrL1, posL1, chrL2, posL2):
    ''' Choose one reprentation of jxn to merge on'''
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


def normalize_jxns(chrom1, chrom2, pos1, pos2, strand1, strand2, repleft, repright, order, args):
    flipstr = string.maketrans("-+", "+-")
    if order == 2:
        if strand1 == "-":
            new_pos1 = chrom1 + ":" + str(pos1 - int(repright)) + ":" + strand1.translate(flipstr)
        else:
            new_pos1 = chrom1 + ":" + str(pos1 + int(repright)) + ":" + strand1.translate(flipstr)
        if strand2 == "-":
            new_pos2 = chrom2 + ":" + str(pos2 - int(repright)) + ":" + strand2.translate(flipstr)
        else:
            new_pos2 = chrom2 + ":" + str(pos2 + int(repright)) + ":" + strand2.translate(flipstr)
        newid = new_pos2 + ":" + new_pos1 + ":" + str(repleft) + ":" + str(repright)
    if order == 1 or args.nucleic_type == "RNA":
        new_pos1 = chrom1 + ":" + str(pos1) + ":" + strand1
        new_pos2 = chrom2 + ":" + str(pos2) + ":" + strand2
        newid = new_pos1 + ":" + new_pos2 + ":" + str(repleft) + ":" + str(repright)
    return newid


def bed_to_tree(bed):
    with open(bed, 'r') as f:
        gtree = GenomeIntervalTree.from_bed(fileobj=f)
    return gtree


def subset_bed_func(jxn, bed_file):
    logger = logging.getLogger("STAR-SEQR")
    logger.info('')
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    targets_tree = bed_to_tree(bed_file)
    intersect1 = len(targets_tree[str(chrom1)].search(int(pos1)))
    intersect2 = len(targets_tree[str(chrom2)].search(int(pos2)))
    if intersect1 or intersect2 >= 1:
        return 1
    else:
        return 0


def get_distance(jxn):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if chrom1 == chrom2:
        dist_between = abs(int(pos1) - int(pos2))
    else:
        dist_between = 10**9
    return dist_between


def get_pairs_func(jxn, rawdf):
    '''
    Get paired end read data that supports each jxn from the junction file.
    discordant paired reads are only represented once.
    These have a -1 for jxntype.
    Need to grab both sides of jxn. Each is unique.
    Later eval the reads from the BAM file.
    This is meant to reduce search space and time processing the BAM later.
    These include duplicates.
    This is slow!!!
    '''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1_less = int(pos1) - 100000
    pos1_plus = int(pos1) + 100000
    pos2_less = int(pos2) - 100000
    pos2_plus = int(pos2) + 100000

    forward = rawdf[(rawdf['jxntype'] == -1) &
                    (rawdf['chrom1'] == chrom1) & (rawdf['chrom2'] == chrom2) &
                    (rawdf['pos1'] > pos1_less) & (rawdf['pos1'] < pos1_plus) &
                    (rawdf['pos2'] > pos2_less) & (rawdf['pos2'] < pos2_plus)]
    reverse = rawdf[(rawdf['jxntype'] == -1) &
                    (rawdf['chrom1'] == chrom2) & (rawdf['chrom2'] == chrom1) &
                    (rawdf['pos1'] > pos2_less) & (rawdf['pos1'] < pos2_plus) &
                    (rawdf['pos2'] > pos1_less) & (rawdf['pos2'] < pos1_plus)]
    return (len(forward), len(reverse))


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
    # parse config
    cfgpaths = parse_config(args.config_file)

    # start analysis
    logger.info("***************STAR-SEQR******************")
    logger.info("CMD = " + str(' '.join(sys.argv)))
    logger.info('Starting to work on sample: ' + args.prefix)

    # run star
    fq1_path = os.path.realpath(args.fastq1)
    fq2_path = os.path.realpath(args.fastq2)
    if args.bed_file:
        bed_path = os.path.realpath(args.bed_file)
    # make sample folder
    if not os.path.exists(args.prefix + "_STAR-SEQR"):
        os.makedirs(args.prefix + "_STAR-SEQR")
    os.chdir(args.prefix + "_STAR-SEQR")
    aux.check_file_exists(fq1_path)
    aux.check_file_exists(fq2_path)
    aux.run_star(cfgpaths, fq1_path, fq2_path, args)
    # start a stats file
    stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')

    # import all jxns
    rawdf = import_jxns_pandas(args.prefix + ".Chimeric.out.junction")
    jxns = rawdf[rawdf['jxntype'] >= 0].reset_index()
    logger.info('Ordering junctions')
    jxns['order'] = jxns.apply(lambda x: choose_order(x['chrom1'], x['pos1'], x['chrom2'], x['pos2']), axis=1)
    # jxns.to_csv(path_or_buf="STAR-SEQR_output_jxns.txt", header=True, sep="\t")
    logger.info('Normalizing junctions')
    jxns['name'] = jxns.apply(lambda x: normalize_jxns(x['chrom1'], x['chrom2'], x['pos1'], x['pos2'],
                                                       x['str1'], x['str2'], x['jxnleft'], x['jxnright'],
                                                       x['order'], args), axis=1)
    groups = jxns.groupby(['name', 'order'], as_index=True)['readid'].agg({'reads': 'count'}).reset_index()
    groups_t = groups.pivot(index='name', columns='order', values='reads').reset_index()
    groups_t.columns = ['name', 'left_counts', 'right_counts']
    # groups_t.to_csv(path_or_buf="STAR-SEQR_output_grouped.txt", header=True, sep="\t")

    # Subset Junctions on ROI, distance, and bidirectionality.
    logger.info('Subsetting junctions')
    if args.bed_file:
        groups_t['subset'] = groups_t.apply(lambda x: subset_bed_func(x['name'], bed_path), axis=1)
        groups_t = groups_t[groups_t['subset'] >= 1]

    if not args.bidir or args.nucleic_type == "RNA":
        tfilt = groups_t[(groups_t["left_counts"] >= args.jxn_reads) |
                         (groups_t["right_counts"] >= args.jxn_reads)].sort_values("left_counts", ascending=False)
    else:
        tfilt = groups_t[(groups_t["left_counts"] >= args.jxn_reads) &
                         (groups_t["right_counts"] >= args.jxn_reads)].sort_values("left_counts", ascending=False)
    # print(tfilt.sort_values("right_counts", ascending=False).head())
    if args.bed_file:
        tfilt['subset'] = tfilt.apply(lambda x: subset_bed_func(x['name'], bed_path), axis=1)
        tfilt = tfilt[tfilt['subset'] >= 1]
    tfilt['dist'] = tfilt.apply(lambda x: get_distance(x['name']), axis=1)

    # Get discordant read pair info. This is slow!
    logger.info('Getting pair info')
    tfilt['pairs_for'], tfilt['pairs_rev'] = zip(*tfilt.apply(lambda x: get_pairs_func(x['name'], rawdf), axis=1))

    logger.info('Filtering junctions based on pairs and distance')
    if not args.bidir or args.nucleic_type == "RNA":
        tfilt2 = tfilt[((tfilt["pairs_for"] >= args.span_reads) |
                        (tfilt["pairs_rev"] >= args.span_reads)) &
                       (tfilt['dist'] >= args.dist)]
    else:
        tfilt2 = tfilt[((tfilt["pairs_for"] >= args.span_reads) &
                        (tfilt["pairs_rev"] >= args.span_reads)) &
                       (tfilt['dist'] >= args.dist)]
    # print(tfilt2.sort_values("right_counts", ascending=False).head())
    tfilt2.to_csv(path_or_buf="STAR-SEQR_candidates.txt", header=True, sep="\t", mode='w',
                  index=False, columns=["name", "left_counts", "right_counts", "pairs_for",
                                        "pairs_rev", "subset", "dist"])

    # Process candidates to see if they pass when duplicates are considered.
    logger.info('Candidates:' + str(len(tfilt2.index)))
    stats_fh.write('Candidates' + '\t' + str(len(tfilt2.index)) + '\n')
    if len(tfilt2.index) >= 1:
        # mark duplicates
        aux.markdups(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.mrkdup.bam", cfgpaths)
        # assemble unique reads
        assemdict = aux.run_assembly_parallel(tfilt2, args.prefix + ".Chimeric.out.mrkdup.bam", cfgpaths, args)
        logger.info("Finished parallel assembly.")
        assemdf = pd.DataFrame.from_records(assemdict, index='name')
        finaldf = pd.merge(tfilt2, assemdf, how='inner', left_on="name", right_on="name", left_index=False,
                           right_index=True, sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
        # finaldf.set_index('name', inplace=True)
        finaldf['jxn_first_unique'] = finaldf["jxnleft_unique_for_first"] + finaldf["jxnleft_unique_rev_first"] + \
            finaldf["jxnright_unique_for_first"] + finaldf["jxnright_unique_rev_first"]
        finaldf['jxn_second_unique'] = finaldf["jxnleft_unique_for_second"] + finaldf["jxnleft_unique_rev_second"] + \
            finaldf["jxnright_unique_for_second"] + finaldf["jxnright_unique_rev_second"]
        finaldf['jxn_first_all'] = finaldf["jxnleft_all_for_first"] + finaldf["jxnleft_all_rev_first"] + \
            finaldf["jxnright_all_for_first"] + finaldf["jxnright_all_rev_first"]
        finaldf['jxn_second_all'] = finaldf["jxnleft_all_for_second"] + finaldf["jxnleft_all_rev_second"] + \
            finaldf["jxnright_all_for_second"] + finaldf["jxnright_all_rev_second"]
        finaldf['ann'] = ann.get_gene_info(cfgpaths['ref_file'], finaldf)
        if not args.bidir or args.nucleic_type == "RNA":
            finaldf = finaldf[(finaldf["spans_disc_unique"] >= args.span_reads) &
                              ((finaldf["jxn_first_unique"] >= args.jxn_reads) |
                               (finaldf["jxn_second_unique"] >= args.jxn_reads))]
        else:
            finaldf = finaldf[(finaldf["spans_disc_unique"] >= args.span_reads) &
                              (finaldf["jxn_first_unique"] >= args.jxn_reads) &
                              (finaldf["jxn_second_unique"] >= args.jxn_reads)]
        # print(finaldf.head)
        finaldf.sort_values(['jxn_first_unique', "spans_disc_unique"], ascending=[False, False], inplace=True)
        outcols = ["ann", "name", "dist", "spans_disc_all", "jxn_first_all", "jxn_second_all",
                   "spans_disc_unique", "jxn_first_unique", "jxn_second_unique",
                   # "jxnleft_unique_for_first", "jxnleft_unique_for_second",
                   # "jxnleft_unique_rev_first", "jxnleft_unique_rev_second",
                   # "jxnright_unique_for_first", "jxnright_unique_for_second",
                   # "jxnright_unique_rev_first", "jxnright_unique_rev_second",
                   "velvet", "primers"]
        if args.spades:
            outcols.append("spades")
            finaldf['primers'] = finaldf.apply(lambda x: primer3.runp3(x['name'], x['spades']), axis=1).apply(lambda x: ",".join(x))
        else:
            finaldf['primers'] = finaldf.apply(lambda x: primer3.runp3(x['name'], x['velvet']), axis=1).apply(lambda x: ",".join(x))
        finaldf.to_csv(path_or_buf=args.prefix + "_STAR-SEQR_breakpoints.txt", header=True, sep="\t",
                       columns=outcols, mode='w', index=False)

        # Make bedpe and VCF
        sv2bedpe.process(finaldf, args)
        logger.info('Breakpoints_identified:' + str(len(finaldf.index)))
        stats_fh.write('Breakpoints' + '\t' + str(len(finaldf.index)) + '\n')
    else:
        logger.info("No candidate junctions identified.")

    # Finish
    end = time.time()
    elapsed = end - start
    logger.info("Program took  %g seconds" % (elapsed))
    # ++ confirm assembled contig is at breakpoint


if __name__ == "__main__":
    sys.exit(main())


# class JXN (object):

#     """object to hold *Chimeric.out.junction from STAR"""

#     def __init__(self, jxnList=[]):
#         if len(jxnList) > 0:
#             self.chr1 = jxnList[0]
#             self.pos1 = int(jxnList[1])
#             self.strand1 = jxnList[2]
#             self.chr2 = jxnList[3]
#             self.pos2 = int(jxnList[4])
#             self.strand2 = jxnList[5]
#             self.jxntype = jxnList[6]
#             self.left_rep = int(jxnList[7])
#             self.right_rep = int(jxnList[8])
#             self.rname = jxnList[9]
#             self.base1_plus = int(jxnList[10])
#             self.cigar1_plus = jxnList[11]
#             self.base2_plus = int(jxnList[12])
#             self.cigar2_plus = jxnList[13]
#             self.valid = 1
#         else:
#             self.valid = 0
#             self.query = 'null'

    # 1: chromosome of the donor
    # 2: First base of the intron of the donor (1-based) (?-1 based)
    # 3: strand of the donor
    # 4: chromosome of the acceptor
    # 5: First base of the intron of the acceptor (1-based) (?-1 based)
    # 6: strand of the acceptor
    # 7: junction type: -1=encompassing junction (between the mates),
    #                                                1=GT/AG, 2=CT/AC
    # 8: repeat length to the left of the junction
    # 9: repeat length to the right of the junction
    #       Columns 10-14 describe the alignments of the two chimeric segments,
    #       it is SAM like. Alignments are given with respect to the (+) strand
    # 10: read name
    # 11: First base of the First segment (on the + strand)
    # 12: CIGAR of the First segment
    # 13: First base of the second segment
    # 14: CIGAR of the second segment

    # SAM
    # 163 flag is _2.fastq.gz read and is full length seq/qual
    # 83 flag is _1.fastq.gz read is reverse/full length of seq/qual
    # 401 flag is chimeric. Seq and Qual is reverse of #163. Cigar is also
    # reverse.

# couple with another method to increase strength

# flags in Test file:
#   11990 65 first, paired, both forward, couples with 129
#   56723 81 first, paired, flip this reads, mate is forward, couples with 161
#    4032 83 proper paired, couples with 163
#   47229 97 first, forward, mate flipped, couples with 145
#    3538 99 proper paired, couples with 147
#   12069 113 first, both flipped, paired, couples with 177
#   11990 129 second, paired, both forward, couples with 65
#   47229 145 second, this is flipped, couples with 97
#    3538 147 proper paied
#   56723 161 second, paired, mate is flipped, this is forward, couples with 81
#    4032 163 proper paied, couples with 83
#   12069 177 second, both flipped, paired, couples with 113
#    1059 321 first, paired, not primary (with 99,147 and 83,163)
#    1187 337 first, paird, reverse strand, not primary (with 83,163)
#    2508 385 second, paired, not primary (with 99,147)
#    2816 401 second, paired, reverse strand, not primary (with, 83,163)
