#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
from six import reraise as raise_
import os
import sys
import re
import logging
import pandas as pd
import numpy as np
import warnings
from collections import OrderedDict
import pysam
import starseqr_utils as su

try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans


logger = logging.getLogger('STAR-SEQR')


def import_starjxns(jxnFile, keep_dups=False):
    ''' read star jxns file and optionally remove dups'''
    logger.info('Importing junctions')
    try:
        df = pd.read_csv(jxnFile, sep="\t", header=None, usecols=range(0, 14), low_memory=False, engine='c')
        df.columns = ['chrom1', 'pos1', 'str1', 'chrom2', 'pos2', 'str2',
                      'jxntype', 'jxnleft', 'jxnright', 'readid',
                      'base1', 'cigar1', 'base2', 'cigar2']
        df['readid'] = df['readid'].astype(str)
        df['pos1'] = df['pos1'].astype(float).astype(int)  # this bypasses some strange numbers
        df['pos2'] = df['pos2'].astype(float).astype(int)
        df['identity'] = df['base1'].astype(str) + ':' + df['cigar1'].astype(str) + ':' + df['base2'].astype(str) + ':' + df['cigar2'].astype(str)
        # df.drop(['base1', 'cigar1', 'base2', 'cigar2'], axis=1, inplace=True)
        if not keep_dups:
            logger.info("Removing duplicate reads")
            return df.drop_duplicates(subset=['identity'], keep='first')
        else:
            logger.info("Using all reads")
            return df
    except Exception as e:
        logger.error("There was a problem reading your STAR *Chimeric.out.junction file")
        logger.error("Exception: " + str(e))
        traceback = sys.exc_info()[2]
        raise_(ValueError, e, traceback)


cigarPattern = '([0-9]+[MIDNSHP])'
cigarSearch = re.compile(cigarPattern)
atomicCigarPattern = '([0-9]+)([MIDNSHP])'
atomicCigarSearch = re.compile(atomicCigarPattern)


def cigar_overhang_matches(cigar1, cigar2):
    matches = 0
    cigar = cigar1
    if "p" in cigar1:
        cigar = cigar2
    if (cigar == "*"):
        return matches
    else:
        cigarsfound = cigarSearch.findall(cigar)
        for opString in cigarsfound:
            cigar_len, cigar_class = atomicCigarSearch.findall(opString)[0]
            cigar_len = int(cigar_len)
            if cigar_class == "M" and cigar_len > matches:
                matches = cigar_len
    return matches


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


def count_jxns(df, nucleic_type="RNA"):
    ''' aggregate jxn reads'''
    grouped_df = df.groupby(['name', 'order'], as_index=True)
    new_df = grouped_df.agg(OrderedDict([('readid', OrderedDict([('reads', lambda col: ','.join(
        col)), ('counts', 'count')])), ('overhang_len', 'max')])).reset_index()
    new_df = new_df.pivot(index='name', columns='order').reset_index()
    if nucleic_type == "DNA":
        new_df.columns = ['name', 'jxnreadsleft', 'jxnreadsright', 'jxnleft', 'jxnright', 'max_overhang']
    elif nucleic_type == "RNA":
        new_df.columns = ['name', 'jxn_reads', 'jxn_counts', 'max_overhang']
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


def exons2seq(fasta_path, lol_exons, jxn, side, fusion_exons='', decorate='', out_dir='support'):
    '''
    Given a fasta and a list of exon coordinates, extract sequence.
    exon boundaries can be annotated with a delimiter of somekind
    header contains: jxn|side|transcripts|pos
    '''
    clean_jxn = su.common.safe_jxn(jxn)
    if out_dir == "support":  # default to support folders
        out_dir = os.path.join('support', clean_jxn)
        out_fa = os.path.join(out_dir, 'transcripts-' + str(side) + ".fa")
    else:
        su.common.make_new_dir(out_dir)
        out_fa = os.path.join(out_dir, 'transcripts-' + str(side) + '-' + clean_jxn + ".fa")

    fa = pysam.Fastafile(fasta_path)
    ofile = open(out_fa, "w")

    all_seq = []
    if sum(1 for x in lol_exons if isinstance(x, list)) == 0:
        ofile.close()
        return  # avoid errors if no exons.

    if fusion_exons:
        if sum(1 for x in fusion_exons if isinstance(x, list)) == 0:
            ofile.close()
            return  # avoid errors if no exons..
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
                        new_trx = str(jxn) + '|' + str(side) + '|' + str(trx) + "--" + (fus_trx) + "|" + str(len(seq_str))
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
                ofile.write(">" + str(jxn) + '|' + str(side) + '|' + ntrx + "\n" + nseq_str + "\n")
            all_seq.append((ntrx, nseq_str))  # need empty if no seq found
    ofile.close()
    return


def mean_from_cols(df, col_regexkey):
    '''
    Function to get mean across a row of specific cols for a df.
    Useful when cols have a comma-sep string of values.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        id_cols = df.filter(regex=col_regexkey, axis=1).columns.tolist()
        res_val = df[id_cols].apply(lambda x: np.nanmean([(float(i) if i else np.nan) for i in ','.join(x.dropna().astype(str)).split(',')]), axis=1)
    return res_val


def minvalcnts_from_cols(df, col_regexkey, minval):
    '''
    Function to get counts greater than a min val across a row of specific cols for a df.
    Useful when cols have a comma-sep string of values.
    '''
    id_cols = df.filter(regex=col_regexkey, axis=1).columns.tolist()
    res_val = df[id_cols].apply(lambda x: np.sum([(float(i) >= minval if i else 0) for i in ','.join(x.dropna().astype(str)).split(',')]), axis=1)
    return res_val


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


def get_svtype_func(jxn):
    '''Used for DNA'''
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
            svtype = "INVERSION"
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


def get_fusion_locations(jxn):
    '''
    STAR junction file reports the first base in the introns, 1-based coordinates
    Returns 0-base coordinates
    '''
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    pos1 = int(pos1) - 1
    pos2 = int(pos2) - 1
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


def write_bedpe(file_in, file_out):
    with open(file_out, 'w') as file_out_fh:
        with open(file_in, 'r') as starOutput:
            for line in starOutput:
                if not line.startswith(('#', 'NAME')):
                    vals = line.split('\t')
                    fusion_name = vals[0]
                    valsL = vals[6].split(':')
                    valsR = vals[7].split(':')
                    chrL = valsL[0]
                    posL = int(valsL[1])
                    strandL = valsL[2]
                    chrR = valsR[0]
                    posR = int(valsR[1])
                    strandR = valsR[2]
                    total_uniqreads = int(vals[1]) + int(vals[2]) + int(vals[3])
                    quant = vals[23]
                    bedpe_line = [chrL, posL, posL + 1,
                                  chrR, posR, posR + 1,
                                  fusion_name, total_uniqreads,
                                  strandL, strandR, quant]

                    bedpe = '\t'.join(list(map(str, bedpe_line)))
                    print(bedpe, file=file_out_fh)


def rna_closeout(prefix, stats_res, breakpoints_fh):
    # Write stats to file
    stats_fh = open(prefix + "_STAR-SEQR.stats", 'w')
    for key, value in stats_res.items():
        stats_fh.write(key + "\t" + str(value) + "\n")
    breakpoints_fh.close()
    brkpt_path = prefix + "_STAR-SEQR_breakpoints.txt"
    bedpe_path = prefix + "_STAR-SEQR_breakpoints.bedpe"
    write_bedpe(brkpt_path, bedpe_path)
