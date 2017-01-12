#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import re
import logging
import pandas as pd
import numpy as np
from collections import OrderedDict
import starseqr_utils as su

try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans


logger = logging.getLogger('STAR-SEQR')


def import_jxns_pandas(jxnFile, args):
    logger = logging.getLogger("STAR-SEQR")
    logger.info('Importing junctions')
    df = pd.read_csv(jxnFile, sep="\t", header=None, usecols=range(0, 14), low_memory=False, engine='c')
    df.columns = ['chrom1', 'pos1', 'str1', 'chrom2', 'pos2', 'str2',
                  'jxntype', 'jxnleft', 'jxnright', 'readid',
                  'base1', 'cigar1', 'base2', 'cigar2']
    df['readid'] = df['readid'].astype(str)
    df['pos1'] = df['pos1'].astype(float).astype(int)  # this bypasses some strange numbers
    df['pos2'] = df['pos2'].astype(float).astype(int)
    df['identity'] = df['base1'].astype(str) + ':' + df['cigar1'].astype(str) + ':' + df['base2'].astype(str) + ':' + df['cigar2'].astype(str)
    df.drop(['base1', 'cigar1', 'base2', 'cigar2'], axis=1, inplace=True)
    if args.keep_dups:
        logger.info("Allowing duplicate reads")
        return df
    else:
        logger.info("Removing duplicate reads")
        return df.drop_duplicates(subset=['identity'], keep='first')


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


def normalize_jxns(chrom1, chrom2, pos1, pos2, strand1, strand2, repleft, repright, order):
    '''Choose one representation for DNA breakpoints'''
    flipstr = maketrans("-+", "+-")
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


def count_jxns(df, args):
    ''' aggregate jxn reads into left and right '''
    grouped_df = df.groupby(['name', 'order'], as_index=True)
    new_df = grouped_df['readid'].agg(OrderedDict([('reads', lambda col: ','.join(col)), ('counts', 'count')]))
    new_df = new_df.reset_index().pivot(index='name', columns='order').reset_index()
    if args.nucleic_type == "DNA":
        new_df.columns = ['name', 'jxnreadsleft', 'jxnreadsright', 'jxnleft', 'jxnright']
    elif args.nucleic_type == "RNA":
        new_df.columns = ['name', 'jxn_reads', 'jxn_counts']
    return new_df


def get_distance(jxn):
    ''' reference distance '''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if chrom1 == chrom2:
        dist_between = abs(int(pos1) - int(pos2))
    else:
        dist_between = np.nan
    return dist_between


def get_pairs_func(jxn, dd):
    '''
    Get paired end read data that supports each jxn from the junction file.
    These have a -1 for jxntype.
    Need to grab both sides of jxn. Each is unique in the .junction file.
    Querying the file can be tedious with the flipping of the junctions.
    Using a dict of pandas dfs to reduce search space and increase performance. 20X increase.
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
    forward = dd[chrom1][(dd[chrom1]['chrom2'] == chrom2) &
                         (dd[chrom1]['pos1'] >= pos1left) & (dd[chrom1]['pos1'] <= pos1right) &
                         (dd[chrom1]['pos2'] >= pos2left) & (dd[chrom1]['pos2'] <= pos2right)]
    reverse = dd[chrom2][(dd[chrom2]['chrom2'] == chrom1) &
                         (dd[chrom2]['pos1'] >= pos2left) & (dd[chrom2]['pos1'] <= pos2right) &
                         (dd[chrom2]['pos2'] >= pos1left) & (dd[chrom2]['pos2'] <= pos1right)]
    npairs = len(forward['readid'].index) + len(reverse['readid'].index)
    reads = ','.join(forward['readid'].tolist()) + ',' + ','.join(reverse['readid'].tolist())

    return (npairs, reads)


def flip_jxn(jxn, gs1):
    '''Flip jxn orientation for RNA Fusion that is inverse according to strand info'''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    chrom1 = str(chrom1)
    chrom2 = str(chrom2)
    pos1 = int(pos1)
    pos2 = int(pos2)
    if str1 != gs1[0]:
        flip = 1
        flipstr = maketrans("-+", "+-")
        if str1 == "-":
            new_pos1 = chrom1 + ":" + str(pos1) + ":" + str1.translate(flipstr)
        else:
            new_pos1 = chrom1 + ":" + str(pos1) + ":" + str1.translate(flipstr)
        if str2 == "-":
            new_pos2 = chrom2 + ":" + str(pos2) + ":" + str2.translate(flipstr)
        else:
            new_pos2 = chrom2 + ":" + str(pos2) + ":" + str2.translate(flipstr)
        newid = new_pos2 + ":" + new_pos1 + ":" + str(repright) + ":" + str(repleft)
    else:
        newid = jxn
        flip = 0
    return (newid, flip)


def exons2seq(fa, lol_exons, jxn, side, fusion_exons='', decorate=''):
    '''
    Given a fasta and a list of exon coordinates, extract sequence.
    exon boundaries can be annotated with a delimiter of somekind
    '''
    clean_jxn = su.common.safe_jxn(jxn)
    jxn_dir = 'support' + '/' + clean_jxn + '/'

    ofile = open(jxn_dir + 'transcripts_' + str(side) + ".fa", "w")
    all_seq = []
    if sum(1 for x in lol_exons if isinstance(x, list)) == 0:
        ofile.close()
        return  # avoid errors if no exons..
    if fusion_exons:
        strand = ''
        fus_strand = ''
        for trx_exons in lol_exons:
            seq = []
            for ex_order, chrom, start, end, strand, trx in trx_exons:  # exon list still ordered in forward
                if fa.references:
                    seq.append(fa.fetch(chrom, start, end))
            if seq:
                seq_str = decorate.join(seq)
                if strand == '-':
                    seq_str = su.common.rc(seq_str)
                for fus_exons in fusion_exons:
                    fus_seq = []
                    for fus_ex_order, fus_chrom, fus_start, fus_end, fus_strand, fus_trx in fus_exons:  # exon list still ordered in forward
                        if fa.references:
                            fus_seq.append(fa.fetch(fus_chrom, fus_start, fus_end))
                    if fus_seq:
                        fus_seq_str = decorate.join(fus_seq)
                        if fus_strand == '-':
                            fus_seq_str = su.common.rc(fus_seq_str)
                        new_trx = str(trx) + "--" + (fus_trx) + "|" + str(len(seq_str))
                        all_seq.append((new_trx, seq_str + fus_seq_str))
                        ofile.write(">" + new_trx + "\n" + seq_str.lower() + fus_seq_str.upper() + "\n")
    else:
        for ntrx_exons in lol_exons:
            nseq = []
            nstrand = ''
            ntrx = ''
            nseq_str = ''
            for nex_order, nchrom, nstart, nend, nstrand, ntrx in ntrx_exons:  # exon list still ordered in forward
                if fa.references:
                    nseq.append(fa.fetch(nchrom, nstart, nend))
            if nseq:
                nseq_str = decorate.join(nseq)
                if nstrand == '-':
                    nseq_str = su.common.rc(nseq_str)
                ofile.write(">" + ntrx + "\n" + nseq_str + "\n")
            all_seq.append((ntrx, nseq_str))  # need empty if no seq found
    ofile.close()
    return pd.Series([[all_seq]])


def get_fusion_class(jxn, txintersection):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    chrom1 = str(chrom1)
    chrom2 = str(chrom2)
    str1 = str(str1)
    str2 = str(str2)

    if len(str(txintersection)) > 3:  # [] counts as two
        return "GENE_INTERNAL"

    if chrom1 != chrom2:
        return "TRANSLOCATION"

    if chrom1 == chrom2:
        if abs(int(pos1) - int(pos2)) < 1000000:
            return "READ_THROUGH"
        elif str1 == str2:
            return "INTERCHROM_INVERTED"
        elif str1 != str2:
            return "INTERCHROM_INTERSTRAND"


def get_minfrag_length(jxn, df):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if str1 == '+':
        left_over = [(int(x) if x else 0) for x in str(df['hangleft_rev_first_seqlen']).split(",")]
        left_jxn = [(int(x) if x else 0) for x in str(df['jxnleft_for_second_seqlen']).split(",")]
    elif str1 == '-':
        left_over = [(int(x) if x else 0) for x in str(df['hangleft_for_first_seqlen']).split(",")]
        left_jxn = [(int(x) if x else 0) for x in str(df['jxnleft_rev_second_seqlen']).split(",")]
    if str2 == '+':
        right_over = [(int(x) if x else 0) for x in str(df['hangright_for_second_seqlen']).split(",")]
        right_jxn = [(int(x) if x else 0) for x in str(df['jxnright_rev_first_seqlen']).split(",")]
    elif str2 == '-':
        right_over = [(int(x) if x else 0) for x in str(df['hangright_rev_second_seqlen']).split(",")]
        right_jxn = [(int(x) if x else 0) for x in str(df['jxnright_for_first_seqlen']).split(",")]
    all_over = left_over + right_over # order matters
    all_jxn = right_jxn + left_jxn # order matters
    all_min = [min(i) for i in list(zip(all_over, all_jxn))] # get the minimum of anchor vs overhang
    all_min20 = np.sum(i > 20 for i in all_min) # how many of the min fragments are > 20
    all_min35 = np.sum(i > 35 for i in all_min) # how many of the min fragments are > 35
    return (all_min20, all_min35)


def get_svtype_func(jxn):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    chrom1 = str(chrom1)
    chrom2 = str(chrom2)
    str1 = str(str1)
    str2 = str(str2)
    if chrom1 != chrom2:
        svtype = "TRANSLOCATION"
    else:
        # STAR notation is same as other tools after strand2 is flipped.
        if str1 == "+" and str2 == "-":
            svtype = "INSERTION"
        elif str1 == "-" and str2 == "+":
            svtype = "INVERSION"
        elif str1 == "+" and str2 == "+":
            svtype = "DELETION"
        elif str1 == "-" and str2 == "-":
            svtype = "DUPLICATION"
    return svtype


def get_sv_locations(jxn):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1 = int(pos1) - int(repright)
    pos2 = int(pos2)
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
    brk1 = str(chrom1) + ":" + str(pos1) + ":" + str(str1)
    brk2 = str(chrom2) + ":" + str(pos2) + ":" + str(str2)
    return (brk1, brk2)
