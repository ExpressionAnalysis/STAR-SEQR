#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import time
import re
import string
from ConfigParser import SafeConfigParser
from argparse import ArgumentParser
import pandas as pd
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import gzip
import starseqr_utils
import starseqr_utils.star_funcs as star
import starseqr_utils.assembly_funcs_dna as assem_dna
import starseqr_utils.assembly_funcs_rna as assem_rna
import starseqr_utils.annotate_sv as ann
import starseqr_utils.sv2bedpe as sv2bedpe
import starseqr_utils.run_primer3 as primer3
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
    parser.add_argument('-d', '--dist', type=int, required=False,
                        default=500,
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
                        default="RNA",
                        help='nucleic acid type',
                        choices=["RNA", "DNA"],
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
    parser.add_argument('-c', '--config_file', type=str, required=False,
                        help='Config file of parameters and paths',
                        metavar="Config file, defaults to " + find_resource("starseqr_config.ini"),
                        default=find_resource("starseqr_config.ini"))
    parser.add_argument('-w', '--workers', type=int, required=False,
                        default=10,
                        help='number of workers',
                        metavar="number of workers")
    parser.add_argument('--spades', '--spades', action='store_true',
                        help='do spades assembly')
    parser.add_argument('--ann_source', '--ann_source', type=str, required=False,
                        default="refgene",
                        help='annotation source',
                        choices=["refgene", "ensgene"],
                        metavar='name of annotaton')
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
    # print(json.dumps(cfgparser._sections['PARAMS'], sort_keys=True, indent=4))
    return cfgparser._sections['PARAMS']


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


def import_jxns_pandas(jxnFile):
    logger = logging.getLogger("STAR-SEQR")
    logger.info('Importing junctions')
    df = pd.read_csv(jxnFile, sep="\t", header=None, usecols=range(0, 10), low_memory=False)
    df.columns = ['chrom1', 'pos1', 'str1', 'chrom2', 'pos2', 'str2', 'jxntype', 'jxnleft', 'jxnright', 'readid']
    df['readid'] = df['readid'].astype(str)
    return df


def choose_order(chrL1, posL1, chrL2, posL2):
    ''' Choose one reprentation of jxn to merge on.
    This is just for DNA breakpoints where bidirectional breakpoint detection occurs
    Also for RNA to flip the name based on gene orientation'''
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
    ''' this is mainly for DNA, but also makes an id for RNA'''
    flipstr = string.maketrans("-+", "+-")
    if order == 2 and not args.nucleic_type == "RNA":
        if strand1 == "-":
            new_pos1 = str(chrom1) + ":" + str(pos1 - int(repright)) + ":" + strand1.translate(flipstr)
        else:
            new_pos1 = str(chrom1) + ":" + str(pos1 + int(repright)) + ":" + strand1.translate(flipstr)
        if strand2 == "-":
            new_pos2 = str(chrom2) + ":" + str(pos2 - int(repright)) + ":" + strand2.translate(flipstr)
        else:
            new_pos2 = str(chrom2) + ":" + str(pos2 + int(repright)) + ":" + strand2.translate(flipstr)
        newid = new_pos2 + ":" + new_pos1 + ":" + str(repleft) + ":" + str(repright)
    elif order == 1 or args.nucleic_type == "RNA":
        new_pos1 = str(chrom1) + ":" + str(pos1) + ":" + strand1
        new_pos2 = str(chrom2) + ":" + str(pos2) + ":" + strand2
        newid = new_pos1 + ":" + new_pos2 + ":" + str(repleft) + ":" + str(repright)
    return newid


def normalize_jxns_rna(chrom1, chrom2, pos1, pos2, strand1, strand2, repleft, repright, order, args):
    ''' this is mainly for DNA, but also makes an id for RNA'''
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


def bed_to_tree(bed):
    with open(bed, 'r') as f:
        gtree = GenomeIntervalTree.from_bed(fileobj=f)
    return gtree


def subset_bed_func(jxn, bed_file):
    logger = logging.getLogger("STAR-SEQR")
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    targets_tree = bed_to_tree(bed_file)
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
        dist_between = 10**9  # max value to use between chromosomes.
    return dist_between


def get_pairs_func(jxn):
    '''
    Get paired end read data that supports each jxn from the junction file.
    discordant paired reads are only represented once.
    These have a -1 for jxntype.
    Need to grab both sides of jxn. Each is unique in the .junction file
    Later eval the reads from the BAM file.
    This is meant to reduce search space and time processing the BAM later.
    These include duplicates.
    This is slow!!! -- consider cython or numba
    '''
    start = time.time()
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if str1 == "+":
        pos1left = int(pos1) - 100000
        pos1right = int(pos1) + int(repright) - 1
    elif str1 == "-":
        pos1left = int(pos1) - int(repright) - 1
        pos1right = int(pos1) + 100000
    if str2 == "+":
        pos2left = int(pos2) - int(repright) - 1
        pos2right = int(pos2) + 100000
    elif str2 == "-":
        pos2left = int(pos2) - 100000
        pos2right = int(pos2) + int(repright) - 1

    forward = rawdf[(rawdf['jxntype'] == -1) &
                    (rawdf['chrom1'] == chrom1) & (rawdf['chrom2'] == chrom2) &
                    (rawdf['pos1'] > pos1left) & (rawdf['pos1'] < pos1right) &.06
                    (rawdf['pos2'] > pos2left) & (rawdf['pos2'] < pos2right)]
    reverse = rawdf[(rawdf['jxntype'] == -1) &
                    (rawdf['chrom1'] == chrom2) & (rawdf['chrom2'] == chrom1) &
                    (rawdf['pos1'] > pos2left) & (rawdf['pos1'] < pos2right) &
                    (rawdf['pos2'] > pos1left) & (rawdf['pos2'] < pos1right)]
    # for_reads = ','.join(forward['readid'].tolist())
    # rev_reads = ','.join(reverse['readid'].tolist())
    print("Span support took  %g seconds" % (time.time() - start))
    return len(forward.index) # + len(reverse.index)


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
    star.check_file_exists(fq1_path)
    star.check_file_exists(fq2_path)
    star.run_star(cfgpaths, fq1_path, fq2_path, args)
    # start output files
    stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')
    breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
    breakpoint_cols = ["ann", "name", "dist", "spans_disc_all", "jxn_first_all", "jxn_second_all",
                       "spans_disc_unique", "jxn_first_unique", "jxn_second_unique",
                       "velvet", "primers"]
    if args.spades:
        breakpoint_cols.append("spades")
    print(*breakpoint_cols, sep='\t', file=breakpoints_fh)

    # import all jxns
    global rawdf
    rawdf = import_jxns_pandas(args.prefix + ".Chimeric.out.junction")
    jxns = rawdf[rawdf['jxntype'] >= 0].reset_index()  # junctions can be either 0, 1, 2
    if args.nucleic_type == "RNA":
        # stats dict
        stats_res = {'Total_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}
        jxns['order'] = 1
        jxns['name'] = jxns.apply(lambda x: normalize_jxns(x['chrom1'], x['chrom2'], x['pos1'], x['pos2'],
                                                           x['str1'], x['str2'], x['jxnleft'], x['jxnright'],
                                                           x['order'], args), axis=1)
        jxns.to_csv(path_or_buf="STAR-SEQR_output_FC.txt", header=True, sep="\t")
        jxn_summary = jxns.groupby(['name', 'order'], as_index=True)['readid'].agg(lambda col: ','.join(col)).reset_index()
        jxn_summary['left_counts'] = jxn_summary['readid'].str.split(',').apply(len)
        stats_res['Total_Breakpoints'] = len(jxn_summary.index)
        logger.info('Total Breakpoints:' + str(stats_res['Total_Breakpoints']))
        if args.bed_file:
            jxn_summary['subset'] = jxn_summary.apply(lambda x: subset_bed_func(x['name'], bed_path), axis=1)
            jxn_summary = jxn_summary[jxn_summary['subset'] >= 1]
        jxn_filt = jxn_summary[(jxn_summary["left_counts"] >= args.jxn_reads)].sort_values("left_counts", ascending=False)

        if len(jxn_filt.index) >= 1:
            # ensembl = pyensembl.EnsemblRelease(release=75)
            # jxn_filt['geneinfo'] = jxn_filt.apply(lambda x: get_genes(x['name'], ensembl), axis=1)
            global gtree
            if args.ann_source == "refgene":
                refgene = find_resource("refGene.txt.gz")
                kg = gzip.open(refgene)
                gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.REF_GENE)
            elif args.ann_source == "ensgene":
                ensgene = find_resource("ensGene.txt.gz")
                kg = gzip.open(ensgene)
                gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)
            jxn_filt['txinfo'] = jxn_filt.apply(lambda x: ann.get_jxn_info(x['name'], gtree), axis=1)
            jxn_filt['dist'] = jxn_filt.apply(lambda x: get_distance(x['name']), axis=1)
            jxn_filt2 = jxn_filt[(jxn_filt['dist'] >= args.dist)]
            # print(jxn_filt2.sort_values("right_counts", ascending=False).head())
            jxn_filt2.to_csv(path_or_buf="STAR-SEQR_candidates.txt", header=True, sep="\t", mode='w',
                             index=False)
            # jxn_filt2['span_counts'] = jxn_filt2.apply(lambda x: get_pairs_func(x['name']), axis=1)
            # Log Stats
            stats_res['Candidate_Breakpoints'] = len(jxn_filt2.index)
            logger.info('Candidate Breakpoints:' + str(stats_res['Candidate_Breakpoints']))

            # Process candidates to see if they pass when duplicates are considered.
            if len(jxn_filt2.index) >= 1:
                # mark duplicates
                star.markdups(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.mrkdup.bam", cfgpaths)
                # get support from unique reads
                assemdict = assem_rna.run_support_parallel(jxn_filt2, args.prefix + ".Chimeric.out.mrkdup.bam", cfgpaths, args)
                logger.info("Finished aggregating support.")
                assemdf = pd.DataFrame.from_records(assemdict, index='name')
                assemdf.to_csv(path_or_buf="STAR-SEQR_candidates_assembled.txt", header=True, sep="\t", mode='w', index=True)
                # Merge Stats for use
                finaldf = pd.merge(jxn_filt2, assemdf, how='inner', left_on="name", right_on="name", left_index=False,
                                   right_index=True, sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
                # collapse reads
                finaldf['jxn_first_unique'] = finaldf["jxnleft_unique_for_first"] + finaldf["jxnleft_unique_rev_first"] + \
                    finaldf["jxnright_unique_for_first"] + finaldf["jxnright_unique_rev_first"]
                finaldf['jxn_second_unique'] = finaldf["jxnleft_unique_for_second"] + finaldf["jxnleft_unique_rev_second"] + \
                    finaldf["jxnright_unique_for_second"] + finaldf["jxnright_unique_rev_second"]
                finaldf['jxn_first_all'] = finaldf["jxnleft_all_for_first"] + finaldf["jxnleft_all_rev_first"] + \
                    finaldf["jxnright_all_for_first"] + finaldf["jxnright_all_rev_first"]
                finaldf['jxn_second_all'] = finaldf["jxnleft_all_for_second"] + finaldf["jxnleft_all_rev_second"] + \
                    finaldf["jxnright_all_for_second"] + finaldf["jxnright_all_rev_second"]
                # TODO: confirm breakpoint with bwa or bowtie
                # TODO: probabilistic module to assign quality score

                # remove this once a qual is assigned.
                finaldf = finaldf[(finaldf["spans_disc_unique"] >= args.span_reads) &
                                  ((finaldf["jxn_first_unique"] >= args.jxn_reads) |
                                  (finaldf["jxn_second_unique"] >= args.jxn_reads))]
                # Generate Primers
                if args.spades:
                    finaldf['primers'] = finaldf.apply(lambda x: primer3.runp3(x['name'], x['spades']), axis=1).apply(lambda x: ",".join(x))
                else:
                    finaldf['primers'] = finaldf.apply(lambda x: primer3.runp3(x['name'], x['velvet']), axis=1).apply(lambda x: ",".join(x))
                # Get Annotation
                finaldf['ann'] = ann.get_gene_info(finaldf, gtree)
                # Write output
                finaldf.sort_values(['jxn_first_unique', "spans_disc_unique"], ascending=[False, False], inplace=True)
                finaldf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                               columns=breakpoint_cols, mode='w', index=False)
                breakpoints_fh.close()
                # Log Stats
                stats_res['Passing_Breakpoints'] = len(finaldf.index)
                logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))
            else:
                logger.info("No candidate junctions identified.")
        else:
            logger.info("No junctions found.")

    elif args.nucleic_type == "DNA":
        # stats dict
        stats_res = {'Total_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}
        jxns['order'] = jxns.apply(lambda x: choose_order(x['chrom1'], x['pos1'], x['chrom2'], x['pos2']), axis=1)
        logger.info('Normalizing junctions')
        jxns['name'] = jxns.apply(lambda x: normalize_jxns(x['chrom1'], x['chrom2'], x['pos1'], x['pos2'],
                                                           x['str1'], x['str2'], x['jxnleft'], x['jxnright'],
                                                           x['order'], args), axis=1)
        groups = jxns.groupby(['name', 'order'], as_index=True)['readid'].agg({'reads': 'count'}).reset_index()
        jxn_summary = groups.pivot(index='name', columns='order', values='reads').reset_index()
        jxn_summary.columns = ['name', 'left_counts', 'right_counts']
        # jxn_summary.to_csv(path_or_buf="STAR-SEQR_output_grouped.txt", header=True, sep="\t")
        logger.info('Subsetting junctions')
        if args.bed_file:
            jxn_summary['subset'] = jxn_summary.apply(lambda x: subset_bed_func(x['name'], bed_path), axis=1)
            jxn_summary = jxn_summary[jxn_summary['subset'] >= 1]
        jxn_filt = jxn_summary[(jxn_summary["left_counts"] >= args.jxn_reads) &
                               (jxn_summary["right_counts"] >= args.jxn_reads)].sort_values("left_counts", ascending=False)
        stats_res['Total_Breakpoints'] = len(jxn_filt.index)
        logger.info('Total Breakpoints:' + str(stats_res['Total_Breakpoints']))
        # Get stats
        if len(jxn_filt.index) >= 1:
            jxn_filt['dist'] = jxn_filt.apply(lambda x: get_distance(x['name']), axis=1)
            # Get discordant read pair info. This is slow!
            logger.info('Getting pair info')
            jxn_filt['pairs_for_id'], jxn_filt['pairs_rev_id'] = zip(*jxn_filt.apply(lambda x: get_pairs_func(x['name'], rawdf), axis=1))
            jxn_filt['pairs_for'] = jxn_filt['pairs_for_id'].str.split(',').apply(len)
            jxn_filt['pairs_rev'] = jxn_filt['pairs_rev_id'].str.split(',').apply(len)
            logger.info('Filtering junctions based on pairs and distance')
            jxn_filt2 = jxn_filt[((jxn_filt["pairs_for"] >= args.span_reads) &
                                  (jxn_filt["pairs_rev"] >= args.span_reads)) &
                                 (jxn_filt['dist'] >= args.dist)]
            jxn_filt2.to_csv(path_or_buf="STAR-SEQR_candidates.txt", header=True, sep="\t", mode='w',
                             index=False, columns=["name", "left_counts", "right_counts", "pairs_for",
                                                   "pairs_rev", "subset", "dist"])

            # Log Stats
            stats_res['Candidate_Breakpoints'] = len(jxn_filt2.index)
            logger.info('Candidate Breakpoints:' + str(stats_res['Candidate_Breakpoints']))
            # Process candidates to see if they pass when duplicates are considered.
            if len(jxn_filt2.index) >= 1:
                # mark duplicates
                star.markdups(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.mrkdup.bam", cfgpaths)
                # assemble unique reads
                assemdict = assem_dna.run_support_parallel(jxn_filt2, args.prefix + ".Chimeric.out.mrkdup.bam", cfgpaths, args)
                logger.info("Finished parallel assembly.")
                assemdf = pd.DataFrame.from_records(assemdict, index='name')
                # Merge Stats for use
                finaldf = pd.merge(jxn_filt2, assemdf, how='inner', left_on="name", right_on="name", left_index=False,
                                   right_index=True, sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
                finaldf['jxn_first_unique'] = finaldf["jxnleft_unique_for_first"] + finaldf["jxnleft_unique_rev_first"] + \
                    finaldf["jxnright_unique_for_first"] + finaldf["jxnright_unique_rev_first"]
                finaldf['jxn_second_unique'] = finaldf["jxnleft_unique_for_second"] + finaldf["jxnleft_unique_rev_second"] + \
                    finaldf["jxnright_unique_for_second"] + finaldf["jxnright_unique_rev_second"]
                finaldf['jxn_first_all'] = finaldf["jxnleft_all_for_first"] + finaldf["jxnleft_all_rev_first"] + \
                    finaldf["jxnright_all_for_first"] + finaldf["jxnright_all_rev_first"]
                finaldf['jxn_second_all'] = finaldf["jxnleft_all_for_second"] + finaldf["jxnleft_all_rev_second"] + \
                    finaldf["jxnright_all_for_second"] + finaldf["jxnright_all_rev_second"]
                # TODO: confirm breakpoint with bwa or bowtie
                # TODO: probabilist module to assign quality score

                # remove this once a qual is assigned.
                finaldf = finaldf[(finaldf["spans_disc_unique"] >= args.span_reads) &
                                  (finaldf["jxn_first_unique"] >= args.jxn_reads) &
                                  (finaldf["jxn_second_unique"] >= args.jxn_reads)]
                # Generate Primers
                if args.spades:
                    finaldf['primers'] = finaldf.apply(lambda x: primer3.runp3(x['name'], x['spades']), axis=1).apply(lambda x: ",".join(x))
                else:
                    finaldf['primers'] = finaldf.apply(lambda x: primer3.runp3(x['name'], x['velvet']), axis=1).apply(lambda x: ",".join(x))
                # Get Annotation
                if args.ann_source == "refgene":
                    refgene = find_resource("refGene.txt.gz")
                    finaldf['ann'] = ann.get_gene_info(refgene, finaldf, "refgene")
                elif args.ann_source == "ensgene":
                    ensgene = find_resource("ensGene.txt.gz")
                    finaldf['ann'] = ann.get_gene_info(ensgene, finaldf, "ensgene")
                # Write output
                finaldf.sort_values(['jxn_first_unique', "spans_disc_unique"], ascending=[False, False], inplace=True)
                finaldf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                               columns=breakpoint_cols, mode='w', index=False)
                breakpoints_fh.close()
                # Log Stats
                stats_res['Passing_Breakpoints'] = len(finaldf.index)
                logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))
            else:
                logger.info("No candidate junctions identified.")
        else:
            logger.info("No junctions found.")

    # Write stats to file
    for key, value in stats_res.items():
        stats_fh.write(key + "\t" + str(value) + "\n")

    # Make bedpe and VCF
    sv2bedpe.process(args.prefix + "_STAR-SEQR_breakpoints.txt", args)

    # Finish
    logger.info("Program took  %g seconds" % (time.time() - start))


if __name__ == "__main__":
    sys.exit(main())
