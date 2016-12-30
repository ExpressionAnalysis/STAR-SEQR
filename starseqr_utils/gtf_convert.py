#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import subprocess
import logging
import gzip

logger = logging.getLogger("STAR-SEQR")


# modified from https://github.com/chapmanb/cloudbiolinux/blob/master/utils/prepare_tx_gff.py
def file_exists(fname):
    """Check if a file exists and is non-empty."""
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False


def gtf_to_genepred(gtf, out_file):
    if file_exists(out_file):
        logger.warning("Skipping gtf_to_genepred as files already exist!")
        return out_file
    cmd = "gtfToGenePred -allErrors -ignoreGroupsWithoutExons -genePredExt -geneNameAsName2 {gtf} {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file


def genepred_to_UCSCtable(genepred, out_file):
    # header = ["#bin", "name", "chrom", "strand",
    #           "txStart", "txEnd", "cdsStart", "cdsEnd",
    #           "exonCount", "exonStarts", "exonEnds", "score",
    #           "name2", "cdsStartStat", "cdsEndStat",
    #           "exonFrames"]
    if file_exists(out_file):
        logger.warning("Skipping genepred_to_UCSC_table as files already exist!")
        return out_file
    with open(genepred) as in_handle, gzip.open(out_file, "wt") as out_handle:
        counter = -1
        current_item = None
        for l in in_handle:
            item = l.split("\t")[0]
            if current_item != item:
                current_item = item
                counter = counter + 1
            out_handle.write("\t".join([str(counter), l]))
    return out_file
