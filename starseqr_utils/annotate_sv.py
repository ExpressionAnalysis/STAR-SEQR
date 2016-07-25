#!/usr/bin/env python

from __future__ import print_function
import gzip
import re
import time
import logging
from intervaltree_bio import GenomeIntervalTree, UCSCTable


logger = logging.getLogger("STAR-SEQR")
# data source for refgene
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/knownGene.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/knownGene.txt.gz
# can also use ucsc or ensg


def get_gene_info(reftable, svtable, kg_type):
    logger.info('Annotating each breakpoint')
    start = time.time()
    kg_open = gzip.open if reftable.endswith('.gz') else open
    kg = kg_open(reftable)
    if kg_type == "refgene":
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.REF_GENE)
    elif kg_type == "ensgene":
        gtree = GenomeIntervalTree.from_table(fileobj=kg, mode='tx', parser=UCSCTable.ENS_GENE)
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
                Lsymbol = list(resL)[idx].data['name2']
                Ltrx = list(resL)[idx].data['name']
                Lstr = list(resL)[idx].data['strand']
                anndict['AnnL'].append((Lsymbol, Ltrx, Lstr))
        else:
            anndict['AnnL'].append(("NA", "NA", "NA"))
        # From the right
        anndict['AnnR'] = []
        if len(resR) > 0:
            for idx, val in enumerate(resR):
                Rsymbol = list(resR)[idx].data['name2']
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
