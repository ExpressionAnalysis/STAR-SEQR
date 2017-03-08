#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import sys
import gzip
import time
import logging
import pandas as pd
import shutil
import subprocess as sp
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')


def cat_all_fa_files(rootDir, out_file, trx_fa=''):
    logger.info("Joining transcript models for salmon")
    fileSet = set()
    for dir_, _, files in os.walk(rootDir):
        for fileName in files:
            if fileName.startswith("transcripts") & fileName.endswith(".fa"):
                relDir = os.path.relpath(dir_, rootDir)
                relFile = os.path.join(rootDir, relDir, fileName)
                fileSet.add(relFile)
    with open(out_file, 'w') as fa_out_fh:
        for f in fileSet:
            with open(f, 'r') as fd:
                shutil.copyfileobj(fd, fa_out_fh)
        if trx_fa:
            if trx_fa.endswith('.gz'):
                trx_fh = gzip.open(trx_fa, 'rt')
            else:
                trx_fh = open(trx_fa, 'r')
            shutil.copyfileobj(trx_fh, fa_out_fh)
            trx_fh.close()
    return


def create_salmon_index(fasta, outdir, nthreads):
    logger.info("Creating an index for salmon")
    su.common.check_file_exists(fasta)
    if os.path.isfile(os.path.join(outdir, "sa.bin")):
        logger.warn("Skipping salmon index as files already exist!")
        return
    cmd = ['salmon', 'index', '-t', fasta, '-i', outdir, '-p', nthreads]
    cmd_args = list(map(str, cmd))
    logger.info("*Command: " + " ".join(cmd_args))
    try:
        p = sp.Popen(cmd_args, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            pass
            # logger.info(stdout)
        if stderr:
            pass
            # logger.error(stderr)
        if p.returncode != 0:
                logger.error("Error: salmon index failed", exc_info=True)
                sys.exit(1)
    except (OSError) as o:
        logger.error("Exception: " + str(o))
        logger.error("salmon index Failed", exc_info=True)
        sys.exit(1)
    return


def run_salmon_quant(index, read1, read2, library, outdir, nthreads):
    logger.info("Running salmon quant")
    if os.path.isfile(os.path.join(outdir, "quant.sf")):
        logger.warn("Skipping salmon quant as files already exist!")
        return
    cmd = ['salmon', 'quant', '-i', index, '-l', library, '-1', read1, '-2', read2,
           '-o', outdir, '-p', nthreads, '-m', 400, '-w', 200, '-q']
    cmd_args = list(map(str, cmd))
    logger.info("*Command: " + " ".join(cmd_args))
    try:
        p = sp.Popen(cmd_args, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if stdout:
            pass
            # logger.debug(stdout)
        if stderr:
            pass
            # logger.error(stderr)
        if p.returncode != 0:
            logger.error("Error: salmon quant failed", exc_info=True)
            sys.exit(1)
    except (OSError) as o:
        logger.error("Exception: " + str(o))
        logger.error("salmon quant Failed", exc_info=True)
        sys.exit(1)
    return


def read_salmon_quant(in_file, all_trx=False):
    logger.info("Summarizing Salmon quant.sf data")
    df = pd.read_csv(in_file, sep="\t", header=0, low_memory=False, engine='c')
    if all_trx:
        df = df[df['Name'].str.contains("fusion|left|right")]
    df['Jxn'], df['Side'], df['Transcript'] = df['Name'].str.split('|', 2).str
    max_df = df.ix[df.groupby(['Jxn', 'Side'], sort=False)['TPM'].idxmax()][['Jxn', 'TPM', 'Side', 'Transcript']].sort_values('Jxn')
    new_df = max_df.pivot(index='Jxn', columns='Side').reset_index()
    new_df.columns = ['Jxn', 'TPM_Fusion', 'TPM_Left', 'TPM_Right', 'Max_Trx_Fusion', 'Max_Trx_Left', 'Max_Trx_Right']
    return new_df


def wrap_salmon(chim_trx_dir, fq1, fq2, library, nthreads, trx_fa=''):
    start = time.time()
    cat_all_fa_files(chim_trx_dir, "candidate_trx_models.fa", trx_fa)
    create_salmon_index("candidate_trx_models.fa", "salmon_index", nthreads)
    run_salmon_quant("salmon_index", fq1, fq2, library, "salmon_quant", nthreads)
    salmon_df = read_salmon_quant("salmon_quant/quant.sf", all_trx=True)
    logger.info('Salmon took  %g seconds' % (time.time() - start))
    return salmon_df
