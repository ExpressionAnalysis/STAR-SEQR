#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import six
import re
import logging
from collections import defaultdict
import pandas as pd


logger = logging.getLogger("STAR-SEQR")

try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans


def convert_str(data):
    if isinstance(data, str):
        return data
    else:
        data_type = type(data)
        if six.PY3 and data_type == six.binary_type:
            return data.decode('utf-8')
        if data_type in (six.string_types, six.text_type):
            return str(data)
        if data_type == dict:
            data = six.iteritems(data)
            return data_type(map(convert_str, data))
        if data_type in (tuple, list, set):
            return data_type(map(convert_str, data))


def clean_intervals(interval_set):
    newres = []
    for idx, val in enumerate(interval_set):
        newres.append(convert_str(list(interval_set)[idx].data))
    return convert_str(newres)


def get_jxn_genes(jxn, gtree):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    resL = clean_intervals(gtree[six.b(chrom1)].search(int(pos1)))
    resR = clean_intervals(gtree[six.b(chrom2)].search(int(pos2)))

    # From the left
    genesL = set()
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            Lsymbol = resL[idx]['name2']
            genesL.add(Lsymbol)
    else:
        genesL.add("NA")
    # From the right
    genesR = set()
    if len(resR) > 0:
        for idx, val in enumerate(resR):
            Rsymbol = resR[idx]['name2']
            genesR.add(Rsymbol)
    else:
        genesR.add("NA")
    return (list(genesL), list(genesR))


def overlap_exon(interval, pos, side):
    '''Input is a single interval... list(resL)[0]'''
    try:
        ends = list(map(int, filter(None, interval['exonEnds'].split(","))))
        starts = list(map(int, filter(None, interval['exonStarts'].split(","))))
        frames = list(map(int, filter(None, interval['exonFrames'].split(","))))
        for idx, x in enumerate(starts):
            if pos in range(starts[idx], ends[idx] + 1):  # end is not inclusive so add 1
                if interval['strand'] == "+":
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
    chrom = interval['chrom']
    starts = list(map(int, filter(None, interval['exonStarts'].split(","))))
    ends = list(map(int, filter(None, interval['exonEnds'].split(","))))
    strand = interval['strand']
    trx = interval['name']
    assert len(ends) == len(starts)
    if coding:
        cs = int(interval['cdsStart'])
        ce = int(interval['cdsEnd'])
    else:
        cs = starts[0]
        ce = ends[-1]
    if brkpt:
        if strand == '+':
            if brkpt_side == 1:
                ce = brkpt
            else:
                cs = brkpt
        else:
            if brkpt_side == 1:
                cs = brkpt
            else:
                ce = brkpt
    if strand == '+':
        pos_list = list(zip(range(1, len(starts) + 1), starts, ends))
        found_exons = [(i, s, e) for i, s, e in pos_list if float(min(e, ce) - max(s, cs)) / (e - s) > 0]
    else:
        pos_list = list(zip(range(len(starts), 0, -1), starts, ends))
        found_exons = [(i, s, e) for i, s, e in pos_list if float(min(e, ce) - max(s, cs)) / (e - s) > 0]
    trimmed_exons = [(i, chrom, max(s, cs), min(e, ce), strand, trx) for i, s, e in found_exons]
    # returns list of order, chrom, start, stop, strand, transcript
    return trimmed_exons


def get_jxnside_anno(jxn, gtree, side):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if side == 2:
        chrom1 = chrom2
        pos1 = pos2
        flipstr = maketrans("-+", "+-")
        str1 = str2.translate(flipstr)
        repleft = repright  # keep to wiggle later
    # convert coordinates to normalized 1-base
    if str1 == "+" and side == 1:
        pos1 = int(pos1) - 1
    elif str1 == "-" and side == 1:
        pos1 = int(pos1)
    elif str1 == "+" and side == 2:
        pos1 = int(pos1) - 1
    elif str1 == "-" and side == 2:
        pos1 = int(pos1)

    resL = clean_intervals(gtree[six.b(chrom1)].search(int(pos1)))
    ann = defaultdict(list)

    if len(resL) > 0:
        for idx, val in enumerate(resL):
            ann['symbol'].append(resL[idx]['name2'])  # symbol
            ann['transcript'].append(resL[idx]['name'])
            ann['strand'].append(resL[idx]['strand'])
            dist, exon, frame = overlap_exon(resL[idx], pos1, side)
            ann['exon'].append(exon)
            ann['dist'].append(dist)
            ann['frame'].append(frame)
            ann['cdslen'].append(int(resL[idx]['cdsEnd']) - int(resL[idx]['cdsStart']))

    else:
        ann['symbol'].append("NA")
        ann['transcript'].append("NA")
        ann['strand'].append("NA")
        ann['exon'].append("NA")
        ann['dist'].append("NA")
        ann['frame'].append("NA")
        ann['cdslen'].append("NA")

    # sort results by dist to exon boundary and then by cdslen and then by transcript for consistency
    dsdf = pd.DataFrame(ann).sort_values(['dist', 'cdslen', 'transcript'], ascending=[True, False, True]).reset_index()
    dsdf['dist'] = dsdf['dist'].astype(str)
    ds = dsdf.to_dict(orient='list')

    ann_string = ','.join([str(a) + ":" + b + ":" + c + ":" +
                           str(d) + ":" + str(e) + ":" + str(f) + ':' +
                           str(g) for a, b, c, d, e, f, g in zip(ds['symbol'], ds['transcript'], ds['strand'], ds['exon'], ds['dist'], ds['frame'], ds['cdslen'])])
    return [ds['symbol'][0], ann_string, ds['strand'][0], ds['cdslen'][0] ]  # just the first values for some fields


def get_pos_genes(chrom1, pos1, gtree):
    resL = clean_intervals(gtree[six.b(chrom1)].search(int(pos1)))
    genesL = set()
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            Lsymbol = resL[idx]['name2']
            genesL.add(Lsymbol)
    else:
        genesL.add("NA")
    return list(genesL)


def get_all_exons(jxn, gtree, side, exon_type="trx"):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if side == 2:
        chrom1 = chrom2
        pos1 = pos2
        flipstr = maketrans("-+", "+-")
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
    resL = clean_intervals(gtree[six.b(chrom1)].search(int(pos1)))
    # From the left

    ann = defaultdict(list)
    if exon_type == "trx":
        if len(resL) > 0:
            for idx, val in enumerate(resL):
                ann['all_exons'].append(get_exons(resL[idx], coding=False))
        else:
            ann['all_exons'].append("NA")
        return ann['all_exons']

    elif exon_type == "fusion":
        if len(resL) > 0:
            for idx, val in enumerate(resL):
                ann['all_exons'].append(get_exons(resL[idx], coding=False, brkpt=pos1, brkpt_side=side))
        else:
            ann['all_exons'].append("NA")
        return ann['all_exons']
    else:
        raise ValueError("Invalid exon_type. Expected one of [trx, fusion]")


def get_gene_region(chrom1, pos1, gtree, pad=500):
    resL = clean_intervals(gtree[six.b(chrom1)].search(int(pos1)))
    gmin = []
    gmax = []
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            gmin.append(int(resL[idx]['txStart']))
            gmax.append(int(resL[idx]['txEnd']))
    region = str(chrom1) + ':' + str(min(gmin) - pad) + '-' + str(max(gmax) + pad)
    return region
