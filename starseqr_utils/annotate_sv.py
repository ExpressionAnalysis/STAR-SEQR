#!/usr/bin/env python

from __future__ import print_function
import re
import logging
from collections import defaultdict
import string
from operator import itemgetter

logger = logging.getLogger("STAR-SEQR")
# data source for refgene
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/knownGene.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/wgEncodeGencodeBasicV24lift37.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/knownGene.txt.gz
# can also use ucsc or ensg


def get_jxn_genes(jxn, gtree):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    resL = gtree[chrom1].search(int(pos1))
    resR = gtree[chrom2].search(int(pos2))
    # From the left
    genesL = set()
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            Lsymbol = list(resL)[idx].data['name2']
            genesL.add(Lsymbol)
    else:
        genesL.add("NA")
    # From the right
    genesR = set()
    if len(resR) > 0:
        for idx, val in enumerate(resR):
            Rsymbol = list(resR)[idx].data['name2']
            genesR.add(Rsymbol)
    else:
        genesR.add("NA")
    return (list(genesL), list(genesR))


def overlap_exon(interval, pos, side):
    '''Input is a single interval... list(resL)[0]'''
    try:
        ends = map(int, filter(None, interval.data['exonEnds'].split(",")))
        starts = map(int, filter(None, interval.data['exonStarts'].split(",")))
        frames = map(int, filter(None, interval.data['exonFrames'].split(",")))
        for idx, x in enumerate(starts):
            if pos in range(starts[idx], ends[idx] + 1):  # end is not inclusive so add 1
                if interval.data['strand'] == "+":
                    if side == 1:
                        dist = ends[idx] - pos
                    else:
                        dist = pos - starts[idx]
                    exon = idx + 1
                else:
                    if side == 1:
                        dist = pos - starts[idx]
                    else:
                        dist = ends[idx] - pos
                    srev = starts[::-1]  # exons are in reverse on "-" strand
                    exon = srev.index(starts[idx]) + 1
                return (dist, exon, frames[idx])
        return("NA", "NA", "NA")
    except:
        return("NA", "NA", "NA")


def get_exons(interval, coding=True, brkpt=False, brkpt_side=1):
    chrom = str(interval.data['chrom'])
    starts = map(int, filter(None, interval.data['exonStarts'].split(",")))
    ends = map(int, filter(None, interval.data['exonEnds'].split(",")))
    strand = interval.data['strand']
    trx = interval.data['name']
    assert len(ends) == len(starts)
    if coding:
        cs = int(interval.data['cdsStart'])
        ce = int(interval.data['cdsEnd'])
    else:
        cs = starts[0]
        ce = ends[-1]
    if brkpt:
        if interval.data['strand'] == '+':
            if brkpt_side == 1:
                ce = brkpt
            else:
                cs = brkpt
        else:
            if brkpt_side == 1:
                cs = brkpt
            else:
                ce = brkpt
    if interval.data['strand'] == '+':
        found_exons = filter(lambda (i, s, e): float(min(e, ce) - max(s, cs)) /
                             (e - s) > 0, zip(range(1, len(starts) + 1), starts, ends))
    else:
        found_exons = filter(lambda (i, s, e): float(min(e, ce) - max(s, cs)) /
                             (e - s) > 0, zip(range(len(starts), 0, -1), starts, ends))
    trimmed_exons = map(lambda (i, s, e): [i, chrom, max(s, cs), min(e, ce), strand, trx], found_exons)
    # returns list of order, chrom, start, stop, strand, transcript
    return trimmed_exons


def get_jxnside_anno(jxn, gtree, side, only_trx=False):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if side == 2:
        chrom1 = chrom2
        pos1 = pos2
        flipstr = string.maketrans("-+", "+-")
        str1 = str2.translate(flipstr)
        repleft = repright  # keep to wiggle later
    if str1 == "+" and side == 1:
        pos1 = int(pos1) - 1
    elif str1 == "-" and side == 1:
        pos1 = int(pos1)
    elif str1 == "+" and side == 2:
        pos1 = int(pos1) - 1
    elif str1 == "-" and side == 2:
        pos1 = int(pos1)
    resL = gtree[chrom1].search(int(pos1))
    # From the left

    ann = defaultdict(list)
    if only_trx:
        if len(resL) > 0:
            for idx, val in enumerate(resL):
                ann['all_exons'].append(get_exons(list(resL)[idx], coding=False))
        else:
            ann['all_exons'].append([])
        return ann['all_exons']

    else:
        if len(resL) > 0:
            for idx, val in enumerate(resL):
                ann['symbol'].append(list(resL)[idx].data['name2'])  # symbol
                ann['transcript'].append(list(resL)[idx].data['name'])
                ann['strand'].append(list(resL)[idx].data['strand'])
                dist, exon, frame = overlap_exon(list(resL)[idx], pos1, side)
                ann['exon'].append(exon)
                ann['dist'].append(dist)
                ann['frame'].append(frame)
                ann['cdslen'].append(int(list(resL)[idx].data['cdsEnd']) - int(list(resL)[idx].data['cdsStart']))
                ann['all_exons'].append(get_exons(list(resL)[idx], coding=False, brkpt=pos1, brkpt_side=side))

        else:
            ann['symbol'].append("NA")
            ann['transcript'].append("NA")
            ann['strand'].append("NA")
            ann['exon'].append("NA")
            ann['dist'].append("NA")
            ann['frame'].append("NA")
            ann['cdslen'].append("NA")
            ann['all_exons'].append("NA")

        # sort results by dist to exon boundary
        ds = {}
        for k in ann.keys():
            _, ds[k] = (list(t) for t in zip(*sorted(zip(ann['dist'], ann[k]), key=itemgetter(0))))

        ann_string = ','.join([str(a) + ":" + b + ":" + c + ":" +
                               str(d) + ":" + str(e) + ":" + str(f) + ':' +
                               str(g) for a, b, c, d, e, f, g in zip(ds['symbol'], ds['transcript'], ds['strand'], ds['exon'], ds['dist'], ds['frame'], ds['cdslen'])])
        # print(ds['symbol'][0], ann_string, ds['strand'][0], ds['cdslen'][0])
        return [ds['symbol'][0], ann_string, ds['strand'][0], ds['cdslen'][0], ds['all_exons']]  # just the first values for some fields


def get_pos_genes(chrom1, pos1, gtree):
    resL = gtree[chrom1].search(int(pos1))
    genesL = set()
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            Lsymbol = list(resL)[idx].data['name2']
            genesL.add(Lsymbol)
    else:
        genesL.add("NA")
    return list(genesL)
