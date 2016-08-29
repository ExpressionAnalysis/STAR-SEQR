#!/usr/bin/env python

from __future__ import print_function
import re
import time
import logging
from collections import defaultdict


logger = logging.getLogger("STAR-SEQR")
# data source for refgene
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/knownGene.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/wgEncodeGencodeBasicV24lift37.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/knownGene.txt.gz
# can also use ucsc or ensg

def get_jxn_info_func(jxn, gtree):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    resL = gtree[chrom1].search(int(pos1))
    resR = gtree[chrom2].search(int(pos2))
    genesL = set()
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            Lsymbol = list(resL)[idx].data['name2']
            genesL.add(Lsymbol)
            # Lstrand = list(resL)[0].data['strand']  # Just use the first gene for now.
    else:
        genesL.add("NA")
    # From the right
    genesR = set()
    if len(resR) > 0:
        for idx, val in enumerate(resR):
            Rsymbol = list(resR)[idx].data['name2']
            genesR.add(Rsymbol)
            # Rstrand = list(resL)[0].data['strand'] # Just use the first gene for now.
    else:
        genesR.add("NA")
    union = list(genesL.union(genesR)) # takes union of genes
    return union


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


def get_gene_info(svtable, gtree):
    logger.info('Annotating each breakpoint')
    start = time.time()
    for index, row in svtable.iterrows():
        chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', row['name'])
        resL = gtree[chrom1].search(int(pos1))
        resR = gtree[chrom2].search(int(pos2))
        # get all possible genes, transcripts that overlap
        # From the left
        anndict = {}
        anndict['AnnL'] = []
        if len(resL) > 0:
            for idx, val in enumerate(resL):
                Lsymbol = list(resL)[idx].data['name2'] # symbol
                Ltrx = list(resL)[idx].data['name']
                Lstr = list(resL)[idx].data['strand']
                anndict['AnnL'].append((Lsymbol, Ltrx, Lstr))
        else:
            anndict['AnnL'].append(("NA", "NA", "NA"))
        # From the right
        anndict['AnnR'] = []
        if len(resR) > 0:
            for idx, val in enumerate(resR):
                Rsymbol = list(resR)[idx].data['name2'] # symbol
                Rtrx = list(resR)[idx].data['name']
                Rstr = list(resR)[idx].data['strand']
                anndict['AnnR'].append((Rsymbol, Rtrx, Rstr))
        else:
            anndict['AnnR'].append(("NA", "NA", "NA"))
        # write to df
        svtable.loc[index, "ann"] = anndict['AnnL'][0][0] + "--" + anndict['AnnR'][0][0]
    end = time.time()
    elapsed = end - start
    logger.info("Annotation took  %g seconds" % (elapsed))
    return(svtable['ann'])


def find_exon(interval, pos, side):
    '''Input is a single interval... list(resL)[0]'''
    try:
        ends = map(int, filter(None, interval.data['exonEnds'].split(",")))
        starts = map(int, filter(None, interval.data['exonStarts'].split(",")))
        frames = map(int, filter(None, interval.data['exonFrames'].split(",")))
        # print(pos)
        # print(*starts)
        # print(*ends)
        for idx, x in enumerate(starts):
            if pos in range(starts[idx], ends[idx] + 1): # end is not inclusive so add 1
                if interval.data['strand'] == "+":
                    if side == 1:
                        dist = ends[idx] - pos
                    else:
                        dist = pos - starts[idx]
                else:
                    if side == 1:
                        dist = pos - starts[idx]
                    else:
                        dist = ends[idx] - pos
                return (dist, idx+1, frames[idx])
        return("NA", "NA", "NA")
    except:
        return("NA", "NA", "NA")


def get_jxnside_anno(jxn, gtree, side):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    if side==2:
        chrom1=chrom2
        pos1=pos2
        str1=str2
        repleft=repright
    # validate this
    if str1=="+" and side==1:
        pos1 = int(pos1) - 1
    elif str1=="-" and side==1:
        pos1 = int(pos1) + 1
    elif str1=="+" and side==2:
        pos1 = int(pos1) + 1
    elif str1=="-" and side==2:
        pos1 = int(pos1) - 1
    resL = gtree[chrom1].search(int(pos1))
    # From the left
    ann = defaultdict(list)
    if len(resL) > 0:
        for idx, val in enumerate(resL):
            ann['symbol'].append(list(resL)[idx].data['name2']) # symbol
            ann['transcript'].append(list(resL)[idx].data['name'])
            ann['strand'].append(list(resL)[idx].data['strand'])
            dist, exon, frame = find_exon(list(resL)[idx], pos1, side)
            ann['exon'].append(exon)
            ann['dist'].append(dist)
            ann['frame'].append(frame)
    else:
        ann['symbol'].append("NA")
        ann['transcript'].append("NA")
        ann['strand'].append("NA")
        ann['exon'].append("NA")
        ann['dist'].append("NA")
        ann['frame'].append("NA")
    ann_string = ','.join([str(a) + ":" + b + ":" + c +":" + str(d) + ":" + str(e) + ":" + str(f) for a,b,c,d,e,f in zip(ann['symbol'],ann['transcript'], ann['strand'], ann['exon'], ann['dist'], ann['frame'])])
    return(ann['symbol'][0], ann_string, ann['strand'][0])  # just the first values for some fields