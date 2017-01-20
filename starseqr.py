#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import sys
import time
import argparse
import pandas as pd
import numpy as np
import logging
from collections import OrderedDict
import starseqr_utils as su


def parse_args():
    class FullPaths(argparse.Action):
        """Expand user- and relative-paths"""
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

    usage = " "
    parser = argparse.ArgumentParser(description="STAR-SEQR Parameters:", epilog=usage)
    # create STAR alignment
    group1 = parser.add_argument_group('Do Alignment', '')
    group1.add_argument('-1', '--fastq1', type=str, required=False, action = FullPaths,
                        help='fastq.gz 1')
    group1.add_argument('-2', '--fastq2', type=str, required=False, action = FullPaths,
                        help='fastq.gz 2')
    group1.add_argument('-i', '--star_index', type=str, required=False, action=FullPaths,
                        help='path to STAR index folder')
    group1.add_argument('-m', '--mode', type=int, required=False,
                        default=0,
                        choices=[0, 1],
                        help='STAR alignment sensitivity Mode. 0=Default, 1=More-Sensitive')

    # existing STAR alignment
    group2 = parser.add_argument_group('Use Existing Alignment', '')
    group2.add_argument('-sj', '--star_jxns', type=str, required=False, action = FullPaths,
                        help='chimeric junctions file produce by STAR')
    group2.add_argument('-ss', '--star_sam', type=str, required=False, action = FullPaths,
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
    parser.add_argument('-r', '--fasta', type=str, required=True, action = FullPaths,
                        help='indexed fasta (.fa|.fa.gz)')
    parser.add_argument('-g', '--gtf', type=str, required=True, action = FullPaths,
                        help='gtf file. (.gtf|.gtf.gz)')
    parser.add_argument('-n', '--nucleic_type', type=str, required=False,
                        default="RNA",
                        help='nucleic acid type',
                        choices=["RNA", "DNA"])
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=8,
                        help='Number of threads to use for STAR and STAR-SEQR. 4-8 recommended.')
    parser.add_argument('-b', '--bed_file', type=str, required=False, action = FullPaths,
                        help='Bed file to subset analysis')
    parser.add_argument('--subset_type', type=str, required=False,
                        default="either",
                        help='allow fusions to pass with either one breakend in bed file or require both. Must use -b.',
                        choices=["either", "both"])
    parser.add_argument('-a', '--as_type', type=str, required=False,
                        default="velvet",
                        help='nucleic acid type',
                        choices=["velvet"])
    parser.add_argument('--keep_dups', action='store_true',
                        help='keep read duplicates. Use for PCR data or with molecular tags')
    parser.add_argument('--keep_gene_dups', action='store_true',
                        help='allow internal gene duplications to be considered')
    parser.add_argument('--keep_mito', action='store_true',
                        help='allow RNA fusions to contain at least one breakpoint from Mitochondria')
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
    _, _, df['test_strand'], _ = zip(*df.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1))
    return df


def apply_flip_func(df):
    df['name'], df['flip'] = zip(*df.apply(lambda x: su.core.flip_jxn(x['name'], x['test_strand']), axis=1))
    return df


def wrap_annotate_junctions(df):
    df['ann_format'] = "Symbol:Transcript:Strand:Exon_No:Dist_to_Exon:Frame:CDS_Length"
    left_ann = df.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 1), axis=1)
    df['left_symbol'], df['left_annot'], df['left_strand'], df['left_cdslen'] = zip(*left_ann)
    right_ann = df.apply(lambda x: su.annotate_sv.get_jxnside_anno(x['name'], gtree, 2), axis=1)
    df['right_symbol'], df['right_annot'], df['right_strand'], df['right_cdslen'] = zip(*right_ann)
    df['ann'] = df['left_symbol'] + "--" + df['right_symbol']
    return df


def wrap_jxn_info(df):
    df['left_canonical'] = df['left_annot'].str.split(':', expand=True)[4].apply(lambda x: 'CANONICAL_SPLICING' if x == '0' else 'NON-CANONICAL_SPLICING')
    df['right_canonical'] = df['right_annot'].str.split(':', expand=True)[4].apply(lambda x: 'CANONICAL_SPLICING' if x == '0' else 'NON-CANONICAL_SPLICING')
    df['splice_type'] = np.where((df['left_canonical'] == 'CANONICAL_SPLICING') &
                                 (df['right_canonical'] == 'CANONICAL_SPLICING'), 'CANONICAL_SPLICING', 'NON-CANONICAL_SPLICING')
    # get all genes associated to look for overlap for each read later..
    df['left_all'], df['right_all'] = zip(*df.apply(lambda x: su.annotate_sv.get_jxn_genes(x['name'], gtree), axis=1))
    df['txunion'] = [list(set(a).union(set(b))) for a, b in zip(df.left_all, df.right_all)]
    df['txintersection'] = [list(set(a).intersection(set(b))) for a, b in zip(df.left_all, df.right_all)]
    # Get reference dist
    df['dist'] = df.apply(lambda x: su.core.get_distance(x['name']), axis=1)
    return df


def apply_get_rna_support(args):
    df, in_bam = args
    dict_res = list(df.apply(lambda x: su.support_funcs_rna.get_rna_support(x['name'], x['txunion'], x['supporting_reads'], in_bam, gtree), axis=1))
    newdf = pd.DataFrame.from_records(dict_res, index='name')
    return newdf  # not passed back in same df


def wrap_exons2seq(args):
    df, fa_object = args
    df['left_trx_exons'] = df.apply(lambda x: su.annotate_sv.get_all_exons(x['name'], gtree, 1, exon_type="trx"), axis=1)
    df['right_trx_exons'] = df.apply(lambda x: su.annotate_sv.get_all_exons(x['name'], gtree, 2, exon_type="trx"), axis=1)
    df['left_fusion_exons'] = df.apply(lambda x: su.annotate_sv.get_all_exons(x['name'], gtree, 1, exon_type="fusion"), axis=1)
    df['right_fusion_exons'] = df.apply(lambda x: su.annotate_sv.get_all_exons(x['name'], gtree, 2, exon_type="fusion"), axis=1)
    df.apply(lambda x: su.core.exons2seq(fa_object, x['left_trx_exons'], x['name'], "left"), axis=1)
    df.apply(lambda x: su.core.exons2seq(fa_object, x['right_trx_exons'], x['name'], "right"), axis=1)
    df.apply(lambda x: su.core.exons2seq(fa_object, x['left_fusion_exons'], x['name'], "all_fusions", x['right_fusion_exons']), axis=1)
    df['write_seq'] = "Finished"
    return df['write_seq'] # sequences are written to fasta not passed, but need to pass something


def apply_primers_func(df):
    df['primers'] = df.apply(lambda x: su.run_primer3.wrap_runp3(x['name'], x['assembly_cross_fusions']), axis=1).apply(lambda x: ",".join(x))
    return df


def apply_get_cross_homology(df):
    df['span_homology_score'], df['jxn_homology_score'] = zip(*df.apply(lambda x: su.cross_homology.get_cross_homology(x['name']), axis=1))
    return df


def apply_get_diversity(df):
    df['overhang_diversity'] = df.apply(lambda x: su.overhang_diversity.get_diversity(x['name']), axis=1)
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


def rna_closeout(prefix, stats_res, breakpoints_fh):
    # Write stats to file
    stats_fh = open(prefix + "_STAR-SEQR.stats", 'w')
    for key, value in stats_res.items():
        stats_fh.write(key + "\t" + str(value) + "\n")
    breakpoints_fh.close()
    brkpt_path = prefix + "_STAR-SEQR_breakpoints.txt"
    bedpe_path = prefix + "_STAR-SEQR_breakpoints.bedpe"
    su.core.write_bedpe(brkpt_path, bedpe_path)


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
    logger.info("STAR-SEQR_version = " + str(su.__version__))
    logger.info('Starting to work on sample: ' + args.prefix)

    # check dependent software can be found
    depend_tools = ['samtools', 'bamfilternames', 'velveth', 'gtfToGenePred']
    for tool in depend_tools:
        if not su.common.which(tool):
            logger.error(tool + "exe not found on path! Quitting.")
            sys.exit(1)

    # check files exist and get abs paths. Necessary since we switch folders.
    depend_paths = [args.fasta, args.gtf, args.fastq1, args.fastq2,
                    args.bed_file, args.star_jxns, args.star_sam]
    for f_item in depend_paths:
        if f_item:
            su.common.check_file_exists(f_item)
            logger.info("Found input: %s", f_item)

    # make sample folder
    if not os.path.exists(args.prefix + "_STAR-SEQR"):
        os.makedirs(args.prefix + "_STAR-SEQR")
    os.chdir(args.prefix + "_STAR-SEQR")

    # Do alignment if fastqs
    if args.fastq1:
        su.star_funcs.run_star(args.fastq1, args.fastq2, args)
    elif args.star_jxns:  # symlink existing files into folder if provided
        su.common.force_symlink(args.star_jxns, args.prefix + ".Chimeric.out.junction")
        su.common.force_symlink(args.star_sam, args.prefix + ".Chimeric.out.sam")

    # import all jxns
    rawdf = su.core.import_starjxns(args.prefix + ".Chimeric.out.junction", args.keep_dups)

    # get jxns only
    jxns = rawdf[rawdf['jxntype'] >= 0].reset_index()  # junctions can be either 0, 1, 2

    # Prepare Annotation
    global gtree # necessary to make global for multiprocessing at the moment
    gtree = su.gtf_convert.gtf2tree(args.gtf)

    if args.nucleic_type == "RNA":
        # start output files
        breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
        breakpoint_cols = ["ann", "span_first", "jxn_left", "jxn_right",
                               "Fusion_Class", "splice_type", "breakpoint_left", "breakpoint_right",
                               "left_symbol", "right_symbol", "ann_format", "left_annot", "right_annot",
                               "dist", "assembly", "assembly_cross_disp", "primers", "name",
                               "span_homology_score", "jxn_homology_score", "overhang_diversity",
                               "minfrag20", "minfrag35", "disposition"]
        breakpoint_header = ["NAME", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT",
                             "FUSION_CLASS", "SPLICE_TYPE", "BRKPT_LEFT", "BRKPT_RIGHT",
                             "LEFT_SYMBOL", "RIGHT_SYMBOL", "ANNOT_FORMAT", "LEFT_ANNOT", "RIGHT_ANNOT",
                             "DISTANCE", "ASSEMBLED_CONTIGS", "ASSEMBLY_CROSS_JXN","PRIMERS", "ID",
                             "SPAN_CROSSHOM_SCORE", "JXN_CROSSHOM_SCORE", "OVERHANG_DIVERSITY",
                             "MINFRAG20", "MINFRAG35", "DISPOSITION"]
        print('\t'.join(map(str, breakpoint_header)), file=breakpoints_fh)

        # stats dict
        stats_res = OrderedDict([('All_Breakpoints', 0), ('Candidate_Breakpoints', 0), ('Passing_Breakpoints', 0)])

        if len(jxns.index) == 0:
            logger.info("No junctions found in the input file")
            rna_closeout(args.prefix, stats_res, breakpoints_fh)
            sys.exit(0)

        # Order, Normalize and Aggregate
        logger.info("Ordering junctions")
        jxns['order'] = 1
        logger.info('Normalizing junctions')
        jxns = su.common.pandas_parallel(jxns, apply_normalize_jxns, args.threads)
        logger.info("Getting gene strand and flipping info as necessary")
        jxns = su.common.pandas_parallel(jxns, apply_jxn_strand, args.threads)
        jxns = su.common.pandas_parallel(jxns, apply_flip_func, args.threads)
        logger.info("Aggregating junctions")
        jxn_summary = su.core.count_jxns(jxns, args.nucleic_type)

        # write stats
        stats_res['All_Breakpoints'] = len(jxn_summary.index)
        logger.info('Total Breakpoints:' + str(stats_res['All_Breakpoints']))

        # Get discordant pairs
        logger.info('Getting pair info')
        # break the rawdf into chromosome specific files for quicker lookups
        dd = {}
        for chrom in set(rawdf['chrom1'].unique()) | set(rawdf['chrom2'].unique()):
            dd[chrom] = rawdf[(rawdf['chrom1'] == chrom) & (rawdf['jxntype'] == -1)]
        jxn_filt = su.common.pandas_parallel(jxn_summary, apply_pairs_func, args.threads, dd)

        # Require at least two reads for processing in order to reduce run time
        logger.info('Filtering junctions')
        before_remove = len(jxn_filt.index)
        jxn_filt = jxn_filt[(jxn_filt["spans"].astype(int) + jxn_filt["jxn_counts"].astype(int)) >= 2]
        logger.info("Number of candidates removed due to total read support less than 2: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # Get Annotation info for each junction
            logger.info('Annotating junctions')
            jxn_filt = su.common.pandas_parallel(jxn_filt, wrap_annotate_junctions, args.threads)

            # determine if junction follows canonical splicing at exon junction
            logger.info('Getting junction info')
            jxn_filt = su.common.pandas_parallel(jxn_filt, wrap_jxn_info, args.threads)

            # subset to ROI using bed file if it exists
            if args.bed_file:
                logger.info('Subsetting junctions using the supplied bed file')
                before_remove = len(jxn_filt.index)
                targets_tree = su.common.bed_to_tree(args.bed_file)
                jxn_filt['subset'] = jxn_filt.apply(lambda x: su.common.subset_bed_func(x['name'], targets_tree, sub_style=args.subset_type), axis=1)
                jxn_filt = jxn_filt[jxn_filt['subset'] >= 1]
                logger.info("Number of candidates removed due to bed filter: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # remove internal gene dups unless otherwise requested, also removes overlapping genes
            if not args.keep_gene_dups:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[jxn_filt['txintersection'].astype(str).str.len() < 3]  # [] counts as two
                logger.info("Number of candidates removed due to internal gene duplication filter: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # remove novel genes unless otherwise requested
            # if not args.keep_novel:
            before_remove = len(jxn_filt.index)
            jxn_filt = jxn_filt[((jxn_filt['left_symbol'] != "NA") & (jxn_filt['right_symbol'] != "NA"))] # No Gene Symbol
            jxn_filt = jxn_filt[((jxn_filt['left_annot'].str.split(':', expand=True)[3] != "NA") & # No Exon number
                                (jxn_filt['right_annot'].str.split(':', expand=True)[3] != "NA"))]
            logger.info("Number of candidates removed due to novel gene filter: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # remove mitochondria unless otherwise requested
            if not args.keep_mito:
                before_remove = len(jxn_filt.index)
                jxn_filt = jxn_filt[~jxn_filt['name'].str.contains("chrM|M")]
                logger.info("Number of candidates removed due to Mitochondria filter: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # remove non-canonical jxns with less than 3 reads here to reduce run time
            jxn_filt = jxn_filt[~((jxn_filt['splice_type'] == "NON-CANONICAL_SPLICING") &
                                   ((jxn_filt["jxn_counts"].astype(int)) < 3))]

        if len(jxn_filt.index) >= 1:
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
            su.common.pandas_parallel(finaldf, wrap_exons2seq, args.threads, args.fasta)  # writes fasta for each transcript per side in support folder

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
            finaldf['assembly_cross_disp'] = finaldf['assembly_cross_fusions'].apply(lambda x: True if len(str(x)) > 1 else False)

            # Generate Primers
            logger.info("Generating primers using indexed fasta")
            finaldf = su.common.pandas_parallel(finaldf, apply_primers_func, args.threads)

            # Get breakpoint locations
            logger.info("Getting normalized breakpoint locations")
            finaldf['breakpoint_left'], finaldf['breakpoint_right'] = zip(*finaldf.apply(lambda x: su.core.get_fusion_locations(x['name']), axis=1))

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
            finaldf['filter_minreads'] = (((finaldf["jxn_first"] >= 2) |  # most robust cases
                                   ((finaldf["span_first"] >= 1) & (finaldf["jxn_left"] >= 1)) |  # read diversity
                                   ((finaldf["span_first"] >= 1) & (finaldf["jxn_right"] >= 1)) |
                                   ((finaldf["jxn_right"] >= 1) & (finaldf["jxn_left"] >= 1))))
            finaldf['filter_minreads'].replace(to_replace=[False], value='minreads', inplace=True, method=None)

            # Hard filter on homology for discordant pairs and jxn. Junctions sequences are usually smaller. Consider a ratio of score to read len?
            finaldf['filter_homology'] = ((finaldf['span_homology_score'] < .40) &
                                          (finaldf['jxn_homology_score'] < .40))
            finaldf['filter_homology'].replace(to_replace=[False], value='homology', inplace=True, method=None)


            # Hard filter on unique overhangs. Requre at least 20% of overhangs to be unique
            finaldf['filter_diversity'] = ((finaldf['overhang_diversity'] >= finaldf['jxn_first'] * .2))
            finaldf['filter_diversity'].replace(to_replace=[False], value='diversity', inplace=True, method=None)

            # Hard filter to require non-canonical splicing events to have greater read support and at least 10% with minfrag35
            noncan_mask = finaldf[finaldf['splice_type'] == "NON-CANONICAL_SPLICING"]
            finaldf['filter_noncanonical'] = ((noncan_mask['jxn_first'] >= 5) &
                                              (noncan_mask['minfrag35'] >= noncan_mask['jxn_first'] * .1))
            finaldf['filter_noncanonical'].replace(to_replace=[False], value='noncanonical_support', inplace=True, method=None)

            # Hard filter to require at least 10% of reads to pass minfrag20 if span reads == 0
            nospan_mask = finaldf[finaldf["span_first"] == 0]
            finaldf['filter_nospanminfrag'] = (nospan_mask['minfrag20'] >= nospan_mask['jxn_first'] * .1)
            finaldf['filter_nospanminfrag'].replace(to_replace=[False], value='nospan_minfrag', inplace=True, method=None)

            # Hard filter to require at least 10% of reads to pass minfrag20 if jxn reads > 10
            highjxn_mask = finaldf[finaldf["jxn_first"] >= 10]
            finaldf['filter_minfrag'] = (highjxn_mask['minfrag20'] >= highjxn_mask['jxn_first'] * .1)
            finaldf['filter_minfrag'].replace(to_replace=[False], value='minfrag', inplace=True, method=None)

            # Get the final disposition of filtering
            finaldf['filter_all'] = finaldf[['filter_minreads', 'filter_homology', 'filter_diversity', 'filter_noncanonical', 'filter_nospanminfrag', 'filter_minfrag']].values.tolist()
            finaldf['disposition'] = finaldf['filter_all'].apply(lambda x: ','.join(x for x in list(map(str, x)) if x not in ['True', 'nan']))
            finaldf['disposition'].replace('', 'PASS', inplace=True)

            # write candidates and info to file
            candid_fh = open(args.prefix + "_STAR-SEQR_candidates.txt", 'w')
            print('\t'.join(map(str, breakpoint_header)), file=candid_fh)
            finaldf.to_csv(path_or_buf=candid_fh, header=False, sep="\t", columns=breakpoint_cols, mode='w', index=False)
            candid_fh.close()

            resultsdf = finaldf[finaldf['disposition'] == 'PASS'].sort_values(['jxn_first', "span_first"], ascending=[False, False])

            # Write passing fusions to file
            resultsdf.to_csv(path_or_buf=breakpoints_fh, header=False, sep="\t",
                             columns=breakpoint_cols, mode='w', index=False)

            # Log Stats
            stats_res['Passing_Breakpoints'] = len(resultsdf.index)
            logger.info('Passing Breakpoints:' + str(stats_res['Passing_Breakpoints']))

            # closeout
            rna_closeout(args.prefix, stats_res, breakpoints_fh)

        else: # there were no candidates to process but still want result files produced
            # closeout
            logger.info('No candidates left to process.')
            rna_closeout(args.prefix, stats_res, breakpoints_fh)


    elif args.nucleic_type == "DNA":
        # start output files
        stats_fh = open(args.prefix + "_STAR-SEQR.stats", 'w')
        breakpoints_fh = open(args.prefix + "_STAR-SEQR_breakpoints.txt", 'w')
        breakpoint_cols = ["ann", "svtype", "breakpoint_left", "breakpoint_right", "dist", "spans_disc", "jxn_first", "jxn_second", "name"]
        breakpoint_header = ["NAME", "SVTYPE", "BRKPT_LEFT", "BRKPT_RIGHT", "DISTANCE", "NREAD_SPANS", "NREAD_JXNLEFT", "NREAD_JXNRIGHT", "ID"]
        print(*breakpoint_header, sep='\t', file=breakpoints_fh)

        # stats dict
        stats_res = {'All_Breakpoints': 0, 'Candidate_Breakpoints': 0, 'Passing_Breakpoints': 0}

        if len(jxns.index) == 0:
            logger.info("No junctions found in the input file")
            su.sv2bedpe.process(jxns, args)
            sys.exit(0)

        # Order, Normalize and Aggregate
        logger.info("Ordering junctions")
        jxns = su.common.pandas_parallel(jxns, apply_choose_order, args.threads)
        logger.info('Normalizing junctions')
        jxns = su.common.pandas_parallel(jxns, apply_normalize_jxns, args.threads)
        logger.info("Aggregating junctions")
        jxn_summary = su.core.count_jxns(jxns, nucleic_type = "DNA")
        # get distance
        jxn_summary['dist'] = jxn_summary.apply(lambda x: su.core.get_distance(x['name']), axis=1)

        # subset to bed if specified
        if args.bed_file:
            start_bed = time.time()
            logger.info('Subsetting junctions using the supplied bed file')
            targets_tree = su.common.bed_to_tree(bed_file)
            jxn_summary['subset'] = jxn_summary.apply(lambda x: su.common.subset_bed_func(x['name'], targets_tree), axis=1)
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
                jxn_filt = jxn_filt[jxn_filt['common'].astype(str).str.len() < 3] # [] counts as two
                logger.info("Number of candidates removed due to internal gene duplication filter: " + str(before_remove - len(jxn_filt.index)))

        if len(jxn_filt.index) >= 1:
            # Get paired discordant spanning reads supporting junctions and filter
            logger.info('Getting pair info')
            dd = {}
            for chrom in set(rawdf['chrom1'].unique()) | set(rawdf['chrom2'].unique()):
                dd[chrom] = rawdf[(rawdf['chrom1'] == chrom) & (rawdf['jxntype'] == -1)]
            jxn_filt = su.common.pandas_parallel(jxn_filt, apply_pairs_func, args.threads, dd)
            logger.info('Filtering junctions based on pairs')
            jxn_filt = jxn_filt[(jxn_filt['spans'] >= args.span_reads)]

            # combine all supporting reads together.
            jxn_filt['supporting_reads'] = jxn_filt['spanreads'] + ',' + jxn_filt['jxnreadsleft'] + ',' + jxn_filt['jxnreadsright']

            # Log Stats
            stats_res['Candidate_Breakpoints'] = len(jxn_filt.index)
            logger.info('Candidate Breakpoints:' + str(stats_res['Candidate_Breakpoints']))

        if len(jxn_filt.index) >= 1:
            # skip the big functions for now
            finaldf = jxn_filt

            finaldf['jxn_first'] = finaldf['jxnleft']
            finaldf['jxn_second'] = finaldf['jxnright']
            finaldf['spans_disc'] = finaldf['spans']

            su.common.sam_2_coord_bam(args.prefix + ".Chimeric.out.sam", args.prefix + ".Chimeric.out.bam", args.threads)
            su.common.check_file_exists(args.prefix + ".Chimeric.out.bam")

            # Get Annotation

            finaldf['ann'] = finaldf['genesleft'].apply(lambda x: str(x[0])) + "--" + finaldf['genesright'].apply(lambda x: str(x[0]))

            # Get Breakpoint type
            finaldf['svtype'] = finaldf.apply(lambda x: su.core.get_svtype_func(x['name']), axis=1)

            # Get breakpoint locations
            finaldf['breakpoint_left'], finaldf['breakpoint_right'] = zip(*finaldf.apply(lambda x: su.core.get_sv_locations(x['name']), axis=1))
            finaldf.to_csv(path_or_buf="STAR-SEQR_candidate_info.txt", header=True, sep="\t", mode='w', index=False)

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
