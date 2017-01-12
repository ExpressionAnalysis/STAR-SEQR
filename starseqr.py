#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import sys
import time
import re
from argparse import ArgumentParser
import pandas as pd
from intervaltree_bio import GenomeIntervalTree, UCSCTable
import gzip
import logging
import numpy as np
import pysam
import starseqr_utils as su


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

    # DNA Parameters
    group3 = parser.add_argument_group('DNA Parameters', '')
    group3.add_argument('-d', '--dist', type=int, required=False,
                        default=100000,
                        help='DNA Only: minimum distance to call junctions')
    group3.add_argument('-j', '--jxn_reads', type=int, required=False,
                        default=2,
                        help='DNA Only: minimum number of junction reads to keep junctions.')
    group3.add_argument('-s', '--span_reads', type=int, required=False,
                        default=2,
                        help='DNA Only: minimum number of spanning discordant read pairs to call junctions.')

    # shared args
    parser.add_argument('-p', '--prefix', type=str, required=True,
                        help='prefix to name files')
    parser.add_argument('-r', '--fasta', type=str, required=True,
                        help='indexed fasta')
    parser.add_argument('-g', '--gtf', type=str, required=True,
                        help='gtf file')
    parser.add_argument('-n', '--nucleic_type', type=str, required=False,
                        default="RNA",
                        help='nucleic acid type',
                        choices=["RNA", "DNA"])
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=8,
                        help='Number of threads to use for STAR and STAR-SEQR. 4-8 recommended.')
    parser.add_argument('-b', '--bed_file', type=str, required=False,
                        help='Bed file to subset analysis')
    parser.add_argument('--subset', type=str, required=False,
                        default="either",
                        help='allow fusions to pass with either one breakend in bed file or require both. Must use -b.',
                        choices=["either", "both"])
    parser.add_argument('-a', '--as_type', type=str, required=False,
                        default="velvet",
                        help='nucleic acid type',
                        choices=["velvet"])
    parser.add_argument('--keep_dups', action='store_true',
                        help='keep read duplicates. Use for PCR data.')
    parser.add_argument('--keep_gene_dups', action='store_true',
                        help='allow internal gene duplications to be considered')
    parser.add_argument('--keep_mito', action='store_true',
                        help='allow RNA fusions to contain at least one breakpoint from Mitochondria')
    # parser.add_argument('--keep_unscaffolded', action='store_true',
    #                    help='allow gene fusions to pass that are found on unscaffolded contigs')
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


def set_log_level_from_verbose(ch, args):
    if not args.verbose:
        ch.setLevel('ERROR')
    elif args.verbose == 1:
        ch.setLevel('WARNING')
    elif args.verbose == 2:
        ch.setLevel('INFO')
    elif args.verbose >= 3:
        ch.setLevel('DEBUG')


def bed_to_tree(bed):
    with open(bed, 'r') as f:
        gtree = GenomeIntervalTree.from_bed(fileobj=f)
    return gtree


def subset_bed_func(jxn, targets_tree, sub_style='either'):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    intersect1 = len(targets_tree[str(chrom1)].search(int(pos1)))
    intersect2 = len(targets_tree[str(chrom2)].search(int(pos2)))
    if sub_style == 'either':
        if intersect1 or intersect2 >= 1:
            return 1
        else:
            return 0
    elif sub_style == 'both':
        if intersect1 or intersect2 >= 1:
            return 1
        else:
            return 0


def apply_choose_order(df):
    df['order'] = df.apply(lambda x: su.core.choose_order(x['chrom1'], x['pos1'], x['chrom2'], x['pos2']), axis=1)
    return df


def apply_normalize_jxns(df):
    df['name'] = df.apply(lambda x: su.core.normalize_jxns(x['chrom1'], x['chrom2'], x['pos1'], x['pos2'],
                                                           x['str1'], x['str2'], x['jxnleft'], x['jxnright'],
                                                           x['order']), axis=1)
    return df


def apply_count_vals(df):
    new_df = df.groupby(['name', 'order'], as_index=True).agg({'readid': {'reads': lambda col: ','.join(
        col), 'counts': 'count'}}).reset_index().pivot(index='name', columns='order').reset_index()
    return new_df


def apply_pairs_func(args):
    df, dd = args
    df['spans'], df['spanreads'] = zip(*df.apply(lambda x: su.core.get_pairs_func(x['name'], dd), axis=1))
    return df


def apply_jxn_strand(df):
    _, _, df['test_strand'], _, _ = zip(*df.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
    return df


def apply_flip_func(df):
    df['name'], df['flip'] = zip(*df.apply(lambda x: su.core.flip_jxn(x['name'], x['test_strand']), axis=1))
    return df


def apply_get_rna_support(args):
    df, in_bam = args
    dict_res = list(df.apply(lambda x: su.support_funcs_rna.get_rna_support(x['name'], x['txunion'], x['supporting_reads'], in_bam, gtree), axis=1))
    newdf = pd.DataFrame.from_records(dict_res, index='name')
    return newdf  # not passed back in same df


def apply_exons2seq(args):
    df, fa_path = args
    fa_object = pysam.Fastafile(fa_path)
    df.apply(lambda x: su.core.exons2seq(fa_object, x['left_trx_exons'], x['name'], "left"), axis=1)
    df.apply(lambda x: su.core.exons2seq(fa_object, x['right_trx_exons'], x['name'], "right"), axis=1)
    df.apply(lambda x: su.core.exons2seq(fa_object, x['left_exons'], x['name'], "all_fusions", x['right_exons']), axis=1)
    return df  # sequences are written to fasta not passed


def apply_primers_func(df):
    df['primers'] = df.apply(lambda x: su.run_primer3.wrap_runp3(x['name'], x['assembly_cross_fusions']), axis=1).apply(lambda x: ",".join(x))
    return df


def apply_get_cross_homology(df):
    df['span_homology_score'], df['jxn_homology_score'] = zip(*df.apply(lambda x: su.cross_homology.get_cross_homology(x['name']), axis=1))
    return df


def apply_get_diversity(df):
    df['overhang_diversity_left'], df['overhang_diversity_right'] = zip(*df.apply(lambda x: su.overhang_diversity.get_diversity(x['name']), axis=1))
    return df


def apply_get_minfrag_length(df):
    df['minfrag20'], df['minfrag35'] = zip(*df.apply(lambda x: su.core.get_minfrag_length(x['name'], x), axis=1))
    return df


def apply_get_assembly_info(args):
    df, as_type = args
    df['assembly'], df['assembly_len'], df['assembly_cross_fusions'] = zip(
        *df.apply(lambda x: su.run_assembly.get_assembly_info(x['name'], as_type), axis=1))
    return df


def apply_get_fusion_class(df):
    df['Fusion_Class'] = df.apply(lambda x: su.core.get_fusion_class(x['name'], x['txintersection']), axis=1)
    return df


def apply_get_annot_db(args):
    df, dbargs = args
    chimerdb3, fuca = dbargs
    df['ChimerDB_ann'] = df.apply(lambda x: su.annotate_db.get_chimerdb(x['name'], chimerdb3), axis=1)
    df['FusionCancer_ann'] = df.apply(lambda x: su.annotate_db.get_fusioncancerdb(x['name'], fuca), axis=1)
    return df


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

    # check dependent software can be found
    if not su.common.which("bamfilternames"):
        logger.error("bamfilternames exe not found on path! Quitting.")
        sys.exit(1)
    if not su.common.which("samtools"):
        logger.error("samtools exe not found on path! Quitting.")
        sys.exit(1)
    if not su.common.which("velveth"):
        logger.error("velveth exe not found on path! Quitting.")
        sys.exit(1)
    if not su.common.which('gtfToGenePred'):
        logger.error("gtfToGenePred not found on path! Quitting.")
        sys.exit(1)

    # check files exist and get abs paths
    if args.fasta:
        fasta_path = os.path.realpath(args.fasta)
        su.common.check_file_exists(fasta_path)
    if args.gtf:
        gtf_path = os.path.realpath(args.gtf)
        su.common.check_file_exists(gtf_path)
    if args.bed_file:
        bed_path = os.path.realpath(args.bed_file)
        su.common.check_file_exists(bed_path)
    if args.fastq1:
        fq1_path = os.path.realpath(args.fastq1)
        fq2_path = os.path.realpath(args.fastq2)
        su.common.check_file_exists(fq1_path)
        su.common.check_file_exists(fq2_path)
    if args.star_jxns:
        starjxns_path = os.path.realpath(args.star_jxns)
        starsam_path = os.path.realpath(args.star_sam)
        su.common.check_file_exists(starjxns_path)
        su.common.check_file_exists(starsam_path)

    # make sample folder
    if not os.path.exists(args.prefix + "_STAR-SEQR"):
        os.makedirs(args.prefix + "_STAR-SEQR")
    os.chdir(args.prefix + "_STAR-SEQR")

    # Do alignment if fastqs
    if args.fastq1:
        su.star_funcs.run_star(fq1_path, fq2_path, args)
    # symlink existing files into folder if provided
    elif args.star_jxns:
        su.common.force_symlink(starjxns_path, args.prefix + ".Chimeric.out.junction")
        su.common.force_symlink(starsam_path, args.prefix + ".Chimeric.out.sam")

    # import all jxns
    rawdf = su.core.import_jxns_pandas(args.prefix + ".Chimeric.out.junction", args)
    jxns = rawdf[rawdf['jxntype'] >= 0].reset_index()  # junctions can be either 0, 1, 2

    if len(jxns.index) == 0:
        logger.info("No junctions found in the input file")
        su.sv2bedpe.process(jxns, args)
        sys.exit(0)

    # Prepare Annotation
    global gtree
    genepred_annot = os.path.splitext(gtf_path)[0] + ".genePred"
    ucsc_annot = os.path.splitext(gtf_path)[0] + ".UCSCTable.gz"
    su.gtf_convert.gtf_to_genepred(gtf_path, genepred_annot)
    su.gtf_convert.genepred_to_UCSCtable(genepred_annot, ucsc_annot)
    kg = gzip.open(ucsc_annot)
    gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)

    if args.nucleic_type == "RNA":
        # start output files
        stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')
        breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
        breakpoint_cols = ["ann", "span_first", "jxn_left", "jxn_right",
                           "Fusion_Class", "splice_type", "breakpoint_left", "breakpoint_right",
                           "left_symbol", "right_symbol", "ann_format", "left_annot", "right_annot",
                           "dist", "assembly", "primers", "name",
                           "span_homology_score", "jxn_homology_score", "overhang_diversity_left", "overhang_diversity_right",
                           "minfrag20", "minfrag35", "PASS"]
        breakpoint_header = ["NAME", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT",
                             "FUSION_CLASS", "SPLICE_TYPE", "BRKPT_LEFT", "BRKPT_RIGHT",
                             "LEFT_SYMBOL", "RIGHT_SYMBOL", "ANNOT_FORMAT", "LEFT_ANNOT", "RIGHT_ANNOT",
                             "DISTANCE", "ASSEMBLY", "PRIMERS", "ID",
                             "SPAN_CROSSHOM_SCORE", "JXN_CROSSHOM_SCORE", "OVERHANG_DIVERSITY_LEFT", "OVERHANG_DIVERSITY_RIGHT",
                             "MINFRAG20", "MINFRAG35", "PASS"]
        print('\t'.join(map(str, breakpoint_header)), file=breakpoints_fh)

        # stats dict
        stats_res = {'All_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}

        # Order, Normalize and Aggregate
        logger.info("Ordering junctions")
        jxns['order'] = 1
        logger.info('Normalizing junctions')
        jxns = su.common.pandas_parallel(jxns, apply_normalize_jxns, args.threads)
        logger.info("Getting gene strand and flipping info as necessary")
        jxns = su.common.pandas_parallel(jxns, apply_jxn_strand, args.threads)
        jxns = su.common.pandas_parallel(jxns, apply_flip_func, args.threads)
        logger.info("Aggregating junctions")
        jxn_summary = su.core.count_jxns(jxns, args)

        # write stats
        stats_res['All_Breakpoints'] = len(jxn_summary.index)
        logger.info('Total Breakpoints:' + str(stats_res['All_Breakpoints']))

        # Get discordant pairs
        logger.info('Getting pair info')
        dd = {}
        for chrom in set(rawdf['chrom1'].unique()) | set(rawdf['chrom2'].unique()):
            dd[chrom] = rawdf[(rawdf['chrom1'] == chrom) & (rawdf['jxntype'] == -1)]
        jxn_filt = su.common.pandas_parallel(jxn_summary, apply_pairs_func, args.threads, dd)

        # hard filter on minimal junction and span reads, require at least two reads...
        logger.info('Filtering junctions')
        before_remove = len(jxn_filt.index)
        jxn_filt = jxn_filt[(jxn_filt["spans"].astype(int) + jxn_filt["jxn_counts"].astype(int)) >= 2]
        logger.info("Number of candidates removed due to total read support less than 2: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # Get dist
            jxn_filt['dist'] = jxn_filt.apply(lambda x: su.core.get_distance(x['name']), axis=1)

            # Get Annotation info for each junction
            jxn_filt['ann_format'] = "Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length"
            jxn_filt['left_symbol'], jxn_filt['left_annot'], jxn_filt['left_strand'], jxn_filt['left_cdslen'], jxn_filt[
                'left_exons'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
            jxn_filt['right_symbol'], jxn_filt['right_annot'], jxn_filt['right_strand'], jxn_filt['right_cdslen'], jxn_filt[
                'right_exons'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 2), axis=1))
            # determine if junction follows canonical splicing at exon junction
            jxn_filt['left_canonical'] = jxn_filt['left_annot'].str.split(':', expand=True)[4].apply(
                lambda x: 'CANONICAL_SPLICING' if x == '0' else 'NON-CANONICAL_SPLICING')
            jxn_filt['right_canonical'] = jxn_filt['right_annot'].str.split(':', expand=True)[4].apply(
                lambda x: 'CANONICAL_SPLICING' if x == '0' else 'NON-CANONICAL_SPLICING')
            jxn_filt['splice_type'] = np.where((jxn_filt['left_canonical'] == 'CANONICAL_SPLICING') &
                                               (jxn_filt['right_canonical'] == 'CANONICAL_SPLICING'),
                                               'CANONICAL_SPLICING', 'NON-CANONICAL_SPLICING')
            # get all genes associated to look for overlap for each read later..
            jxn_filt['left_all'], jxn_filt['right_all'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxn_genes(x['name'], gtree), axis=1))
            jxn_filt['txunion'] = [list(set(a).union(set(b))) for a, b in zip(jxn_filt.left_all, jxn_filt.right_all)]
            jxn_filt['txintersection'] = [list(set(a).intersection(set(b))) for a, b in zip(jxn_filt.left_all, jxn_filt.right_all)]
            jxn_filt['ann'] = jxn_filt['left_symbol'] + "--" + jxn_filt['right_symbol']

            # jxn_filt.to_csv(path_or_buf="transformed.txt", header=True, sep=str('\t'), mode='w', index=False)

            # subset to ROI using bed file if it exists
            if args.bed_file:
                logger.info('Subsetting junctions using the supplied bed file')
                targets_tree = bed_to_tree(bed_path)
                jxn_filt['subset'] = jxn_filt.apply(lambda x: subset_bed_func(x['name'], targets_tree, sub_style=args.subset), axis=1)
                jxn_filt = jxn_filt[jxn_filt['subset'] >= 1]

            # remove internal gene dups unless otherwise requested, also removes overlapping genes
            if not args.keep_gene_dups:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[jxn_filt['txintersection'].astype(str).str.len() < 3]  # [] counts as two
                logger.info("Number of candidates removed due to internal gene duplication filter: " + str(before_remove - len(jxn_filt.index)))

            # remove novel genes unless otherwise requested
            # if not args.keep_novel:
            before_remove = len(jxn_filt.index)
            jxn_filt = jxn_filt[(jxn_filt['left_symbol'] != "NA") & (jxn_filt['right_symbol'] != "NA")]
            jxn_filt = jxn_filt[(jxn_filt['left_annot'].str.split(':', expand=True)[3] != "NA") &
                                (jxn_filt['right_annot'].str.split(':', expand=True)[3] != "NA")]
            logger.info("Number of candidates removed due to novel gene filter: " + str(before_remove - len(jxn_filt.index)))

            # remove mitochondria unless otherwise requested
            if not args.keep_mito:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[~jxn_filt['name'].str.contains("chrM|M")]
                logger.info("Number of candidates removed due to Mitochondria filter: " + str(before_remove - len(jxn_filt.index)))

            # remove non-canonical contigs unless otherwise requested
            # if not args.keep_unscaffolded:
            #     before_remove = len(jxn_filt.index)
            #     jxn_filt = jxn_filt[~jxn_filt['name'].str.contains("Un|random|hap|GL0|NC")]
            #     logger.info("Number of candidates removed due to non-canonical contigs filter: " + str(before_remove - len(jxn_filt.index)))

            # remove non-coding
            # if not args.keep_noncoding:
            #     before_remove = len(jxn_filt.index)
            #     jxn_filt = jxn_filt[(jxn_filt['left_cdslen'] != 0) & (jxn_filt['right_cdslen'] != 0)]
            #     logger.info("Number of candidates removed due to non-coding filter: " + str(before_remove - len(jxn_filt.index)))

            # combine all supporting reads together.
            jxn_filt['supporting_reads'] = jxn_filt['jxn_reads'] + ',' + jxn_filt['spanreads']

            stats_res['Candidate_Breakpoints'] = len(jxn_filt.index)
            logger.info('Candidate Breakpoints:' + str(stats_res['Candidate_Breakpoints']))

        # Process candidates
        if len(jxn_filt.index) >= 1:
            # Convert sam to bam
            su.common.sam_2_coord_bam(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.bam", args.threads)
            su.common.check_file_exists(args.prefix + ".Chimeric.out.bam")

            # Gather unique read support
            logger.info("Getting read support.")
            su.common.make_new_dir('support')
            in_bam = args.prefix + ".Chimeric.out.bam"
            support_df = su.common.pandas_parallel(jxn_filt, apply_get_rna_support, args.threads, in_bam)
            finaldf = pd.merge(jxn_filt, support_df, how='inner', left_on="name", right_on="name", left_index=False,
                               right_index=True, sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
            # collapse read info for brevity but keep here in case useful later on
            finaldf['jxn_left'] = finaldf["jxnleft_for_first"] + finaldf["jxnleft_rev_first"]
            finaldf['jxn_right'] = finaldf["jxnright_for_first"] + finaldf["jxnright_rev_first"]
            finaldf['jxn_first'] = finaldf["jxnleft_for_first"] + finaldf["jxnleft_rev_first"] + \
                finaldf["jxnright_for_first"] + finaldf["jxnright_rev_first"]
            finaldf['jxn_second'] = finaldf["jxnleft_for_second"] + finaldf["jxnleft_rev_second"] + \
                finaldf["jxnright_for_second"] + finaldf["jxnright_rev_second"]
            finaldf['span_first'] = finaldf["spanleft_for_first"] + finaldf["spanleft_rev_first"] + \
                finaldf["spanright_for_first"] + finaldf["spanright_rev_first"]
            finaldf['span_second'] = finaldf["spanleft_for_second"] + finaldf["spanleft_rev_second"] + \
                finaldf["spanright_for_second"] + finaldf["spanright_rev_second"]
            finaldf['spans_disc'] = finaldf['span_first']

            # Get all overlapping transcript seqs into one fasta per side
            finaldf['left_trx_exons'] = finaldf.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1, only_trx=True), axis=1)
            finaldf['right_trx_exons'] = finaldf.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 2, only_trx=True), axis=1)
            finaldf = su.common.pandas_parallel(finaldf, apply_exons2seq, args.threads, fasta_path)  # returns sequences for each transcript per side.

            # get homology mapping scores
            logger.info("Getting read homology mapping scores")
            finaldf = su.common.pandas_parallel(finaldf, apply_get_cross_homology, args.threads)

            # get overhang read diversity
            logger.info("Getting overhang read diversity")
            finaldf = su.common.pandas_parallel(finaldf, apply_get_diversity, args.threads)

            # get min frag lengths for anchor/overhang
            logger.info("Getting min frag length support")
            finaldf = su.common.pandas_parallel(finaldf, apply_get_minfrag_length, args.threads)

            # get assembly seq and confirm breakpoint
            logger.info("doing assembly")
            finaldf = su.common.pandas_parallel(finaldf, apply_get_assembly_info, args.threads, args.as_type)

            # Generate Primers
            logger.info("Generating primers using indexed fasta")
            finaldf = su.common.pandas_parallel(finaldf, apply_primers_func, args.threads)

            # Get breakpoint locations
            logger.info("Getting normalized breakpoint locations")
            finaldf['breakpoint_left'], finaldf['breakpoint_right'] = zip(*finaldf.apply(lambda x: su.core.get_sv_locations(x['name']), axis=1))

            # get homology mapping scores
            logger.info("Getting read homology mapping scores")
            finaldf = su.common.pandas_parallel(finaldf, apply_get_fusion_class, args.threads)

            # Get annotations for known fusions
            logger.info("Getting annotation")
            chimerdb3_path = su.common.find_resource('ChimerDB3.0_ChimerSeq.bedpe')
            fuca_path = su.common.find_resource('FusionCancer.bedpe')
            chimerdb3 = pd.read_csv(chimerdb3_path, sep="\t", header=None, low_memory=False, engine='c')
            chimerdb3.columns = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Fusion_pair', 'Counts',
                                 'strand1', 'strand2', 'Source', 'Cancer', 'Features']
            fuca = pd.read_csv(fuca_path, sep="\t", header=None, low_memory=False, engine='c')
            fuca.columns = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Database_ID', 'Rate',
                            'strand1', 'strand2', 'Methods', 'Cancers', 'Fusion_name']
            finaldf = su.common.pandas_parallel(finaldf, apply_get_annot_db, args.threads, chimerdb3, fuca)

            # TODO: probabilistic module to assign quality score

            # HARD FILTERING - Change this once a probabilistic module is ready.
            # Hard filter on read counts after accounting for transcript info.
            finaldf['filter1'] = (((finaldf["jxn_first"] >= 2) |  # most robust cases
                                   ((finaldf["span_first"] >= 1) & (finaldf["jxn_left"] >= 1)) |  # read diversity
                                   ((finaldf["span_first"] >= 1) & (finaldf["jxn_right"] >= 1)) |
                                   ((finaldf["jxn_right"] >= 1) & (finaldf["jxn_left"] >= 1))))

            # Hard filter on homology for discordant pairs and jxn. Junctions sequences are usually smaller. Consider a ratio of score to read len?
            finaldf['filter2'] = ((finaldf['span_homology_score'] < 40) &
                                  (finaldf['jxn_homology_score'] < 40))

            # Hard filter on unique overhangs. Requre at least 20% of overhangs to be unique
            finaldf['filter3'] = ((finaldf['overhang_diversity_left'] + finaldf['overhang_diversity_right'] >= finaldf['jxn_first'] * .2))

            # Hard filter to require non-canonical splicing events to have greater read support and at least 10% with minfrag35
            finaldf['filter4'] = ((finaldf[finaldf['splice_type'] == "NON-CANONICAL_SPLICING"]['jxn_first'] >= 5) &
                                  (finaldf['minfrag35'] >= finaldf['jxn_first'] * .1))

            # Hard filter to require at least 10% of reads to pass minfrag20 if span reads < 1
            finaldf['filter5'] = (finaldf[finaldf["span_first"] < 1]['minfrag20'] >= finaldf[finaldf["span_first"] < 1]['jxn_first'] * .1)

            # Hard filter to require at least 10% of reads to pass minfrag20 if jxn reads > 10
            finaldf['filter6'] = (finaldf[finaldf["jxn_first"] >= 10]['minfrag20'] >= finaldf[finaldf["jxn_first"] >= 10]['jxn_first'] * .1)

            finaldf['PASS'] = finaldf[['filter1', 'filter2', 'filter3', 'filter4', 'filter5', 'filter6']].all(axis=1)

            # all candidates
            candid_fh = open(args.prefix + "_STAR-SEQR_candidates.txt", 'w')
            candid_cols = ["ann", "span_first", "jxn_left", "jxn_right",
                               "Fusion_Class", "splice_type", "breakpoint_left", "breakpoint_right",
                               "left_symbol", "right_symbol", "ann_format", "left_annot", "right_annot",
                               "dist", "assembly", "primers", "name",
                               "span_homology_score", "jxn_homology_score", "overhang_diversity_left", "overhang_diversity_right",
                               "minfrag20", "minfrag35", "PASS"]
            candid_header = ["NAME", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT",
                                 "FUSION_CLASS", "SPLICE_TYPE", "BRKPT_LEFT", "BRKPT_RIGHT",
                                 "LEFT_SYMBOL", "RIGHT_SYMBOL", "ANNOT_FORMAT", "LEFT_ANNOT", "RIGHT_ANNOT",
                                 "DISTANCE", "ASSEMBLY", "PRIMERS", "ID",
                                 "SPAN_CROSSHOM_SCORE", "JXN_CROSSHOM_SCORE", "OVERHANG_DIVERSITY_LEFT", "OVERHANG_DIVERSITY_RIGHT",
                                 "MINFRAG20", "MINFRAG35", "PASS"]
            print('\t'.join(map(str, candid_header)), file=candid_fh)
            finaldf.to_csv(path_or_buf=candid_fh, header=False, sep="\t", columns=candid_cols, mode='w', index=False)
            candid_fh.close()

            resultsdf = finaldf[finaldf['PASS'] == True].sort_values(['jxn_first', "span_first"], ascending=[False, False])

            # Write output
            resultsdf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                             columns=breakpoint_cols, mode='w', index=False)
            breakpoints_fh.close()

            # Make bedpe and VCF
            su.sv2bedpe.process(resultsdf, args)

            # Log Stats
            stats_res['Passing_Breakpoints'] = len(resultsdf.index)
            logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))

        if len(jxn_filt.index) == 0:
            logger.info("No junctions found after filtering")
            su.sv2bedpe.process(jxn_filt, args)
            sys.exit(0)

    elif args.nucleic_type == "DNA":
        # start output files
        stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')
        breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
        breakpoint_cols = ["ann", "svtype", "breakpoint_left", "breakpoint_right", "dist", "spans_disc", "jxn_first", "jxn_second", "name"]
        breakpoint_header = ["NAME", "SVTYPE", "BRKPT_LEFT", "BRKPT_RIGHT", "DISTANCE", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT", "ID"]
        print(*breakpoint_header, sep='\t', file=breakpoints_fh)

        # stats dict
        stats_res = {'All_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}

        # Order, Normalize and Aggregate
        logger.info("Ordering junctions")
        jxns = su.common.pandas_parallel(jxns, apply_choose_order, args.threads)
        logger.info('Normalizing junctions')
        jxns = su.common.pandas_parallel(jxns, apply_normalize_jxns, args.threads)
        logger.info("Aggregating junctions")
        jxn_summary = su.core.count_jxns(jxns, args)

        # get distance
        jxn_summary['dist'] = jxn_summary.apply(lambda x: su.core.get_distance(x['name']), axis=1)

        # subset to bed if specified
        if args.bed_file:
            start_bed = time.time()
            logger.info('Subsetting junctions using the supplied bed file')
            targets_tree = bed_to_tree(bed_path)
            jxn_summary['subset'] = jxn_summary.apply(lambda x: subset_bed_func(x['name'], targets_tree), axis=1)
            jxn_summary = jxn_summary[jxn_summary['subset'] >= 1]
            logger.info("Time to subset junction from bed took %g seconds" % (time.time() - start_bed))

        # print stats
        stats_res['All_Breakpoints'] = len(jxn_summary.index)
        logger.info('Total Breakpoints:' + str(stats_res['All_Breakpoints']))

        # Filter on distance
        logger.info('Filtering junctions based on distance')
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
            jxn_filt['genesleft'], jxn_filt['genesright'] = zip(*jxn_filt.apply(lambda x: su.annotate_sv.get_jxn_genes(x['name'], gtree), axis=1))
            jxn_filt['common'] = [list(set(a).intersection(set(b))) for a, b in zip(jxn_filt.genesleft, jxn_filt.genesright)]

            # Get gene info and remove internal gene dups
            if not args.keep_gene_dups:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[jxn_filt['common'] == 0]
                logger.info("Number of candidates removed due to internal gene duplication filter: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # Get paired discordant spanning reads supporting junctions and filter
            logger.info('Getting pair info')
            dd = {}
            for chrom in set(rawdf['chrom1'].unique()) | set(rawdf['chrom2'].unique()):
                dd[chrom] = rawdf[(rawdf['chrom1'] == chrom) & (rawdf['jxntype'] == -1)]
            jxn_filt = su.common.pandas_parallel(jxn_summary, apply_pairs_func, args.threads, dd)
            logger.info('Filtering junctions based on pairs')
            jxn_filt = jxn_filt[(jxn_filt['spans'] >= args.span_reads)]

            # combine all supporting reads together.
            jxn_filt['supporting_reads'] = jxn_filt['spanreads'] + ',' + jxn_filt['jxnreadsleft'] + ',' + jxn_filt['jxnreadsright']

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
            su.common.sam_2_coord_bam(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.bam", args.threads)
            su.common.check_file_exists(args.prefix + ".Chimeric.out.bam")

            # Get Annotation
            finaldf['ann'] = finaldf['genesleft'].apply(lambda x: str(x[0])) + "--" + finaldf['genesright'].apply(lambda x: str(x[0]))

            # Get Breakpoint type
            finaldf['svtype'] = finaldf.apply(lambda x: su.core.get_svtype_func(x['name']), axis=1)

            # Get breakpoint locations
            finaldf['breakpoint_left'], finaldf['breakpoint_right'] = zip(*finaldf.apply(lambda x: su.core.get_sv_locations(x['name']), axis=1))

            finaldf.to_csv(path_or_buf="STAR-SEQR_candidate_info.txt", header=True, sep="\t", mode='w', index=False)
            # finaldf = finaldf[(finaldf["spans_disc_unique"] >= args.span_reads) &
            #                   (finaldf["jxn_first_unique"] >= args.jxn_reads) &
            #                   (finaldf["jxn_second_unique"] >= args.jxn_reads)]

            # Write output
            # finaldf['dist'] = finaldf['dist'].str.replace()
            finaldf.sort_values(['jxnleft', "spans"], ascending=[False, False], inplace=True)
            finaldf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                           columns=breakpoint_cols, mode='w', index=False)
            breakpoints_fh.close()

            # Make bedpe and VCF
            su.sv2bedpe.process(finaldf, args)

            # Log Stats
            stats_res['Passing_Breakpoints'] = len(finaldf.index)
            logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))
        else:
            logger.info("No candidate junctions identified.")
            # Make bedpe and VCF, write headers only
            su.sv2bedpe.process(jxn_filt, args)

    # Write stats to file
    for key, value in stats_res.items():
        stats_fh.write(key + "\t" + str(value) + "\n")

    # Finish
    logger.info("Program took  %g seconds" % (time.time() - start))


if __name__ == "__main__":
    sys.exit(main())
