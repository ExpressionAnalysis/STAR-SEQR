#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import gzip
import logging
import time
import errno
import string
import pandas as pd
import numpy as np
import itertools
from itertools import groupby
import multiprocessing as mp
import signal
import pysam
import starseqr_utils

logger = logging.getLogger('STAR-SEQR')


def check_file_exists(path):
    if (os.stat(os.path.realpath(path)).st_size == 0):
        logger.error("Exiting. Cannot find file: " + os.path.realpath(path))
        sys.exit(1)


def is_newer_than(target_path, orig_path):
    if not os.path.isfile(target_path):
        return False
    return (os.stat(target_path).st_mtime >= os.stat(orig_path).st_mtime)


def make_new_dir(newdir):
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass


def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)


def which(program):

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def find_resource(filename):
    packagedir = starseqr_utils.__path__[0]
    dirname = os.path.join(packagedir, 'resources')
    fullpath = os.path.abspath(os.path.join(dirname, filename))
    return fullpath


def pandas_parallel(df, func, nthreads, *opts):
    '''wrapper to run pandas apply func using map. *opts will take any number of optional arguments to be passed
    into the function. Note these will be passed to the function as a tuple and need to be parsed. '''
    start = time.time()
    def init_worker():
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        df_split = np.array_split(df, min(nthreads, len(df.index)))
        pool = mp.Pool(nthreads, init_worker)
        if not opts:
            tuple_opts = df_split
        if len(opts) == 1:
            tuple_opts = itertools.izip(df_split, itertools.repeat(opts[0]))
        elif len(opts) > 1:
            tuple_opts = itertools.izip(df_split, itertools.repeat(opts))
        these_res = pool.map(func, tuple_opts)
        df = pd.concat(these_res)
        pool.close()
        pool.join()
    except KeyboardInterrupt as e:
        logger.error("Error: Keyboard interrupt")
        pool.terminate()
        raise e
    except Exception as e:
        logger.error("Exception: " + str(e))
        pool.terminate()
        raise e
    logger.info("Time to run pandas_parallel on " + str(func.__name__) + " took %g seconds" % (time.time() - start))
    return df


def safe_jxn(jxn):
    ''' make a jxn string safe for folder and filenames'''
    fix_jxn = str(jxn).replace(':', '_')
    fix_jxn = str(fix_jxn).replace('+', 'pos')
    fix_jxn = str(fix_jxn).replace('-', 'neg')
    return fix_jxn


def rc(dna):
    ''' reverse complement '''
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def fasta_iter(fasta_name):
    ''' Given a fasta file, yield tuples of header, sequence'''
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
    for header in faiter:
        header = header.next()[1:].strip() # drop the '>'
        # join all sequence lines to one.
        seq = ''.join(s.strip() for s in faiter.next())
        yield header, seq
    fh.close()


class FastqRead(object):
    """Represents 1 read, or 4 lines of a casava 1.8-style fastq file."""
    def __init__(self, file_obj):
        header = file_obj.next().rstrip()
        assert header.startswith('@')
        self.header = header
        self.sequence = file_obj.next().rstrip()
        self.seq_len = len(self.sequence)
        file_obj.next()
        self.quality = file_obj.next().rstrip()
    def __str__(self):
        return "{0}\n{1}\n+\n{2}\n+\n{3}\n".format(self.header, self.sequence,
                                           self.quality, self.seq_len)


class FastqParser(object):
    """Parses a fastq file into FastqReads."""
    def __init__(self, filename, parse_headers=True):
        if filename.endswith('.gz'):
            self._file = gzip.open(filename, 'rb')
        else:
            self._file = open(filename, 'rU')
        self._line_index = 0
    def __iter__(self):
        return self
    def next(self):
        read = FastqRead(self._file)
        if read:
            self._line_index += 4
            return read
        else:
            self._file.close()
            raise StopIteration()


def sam_2_coord_bam(in_sam, out_bam, nthreads):
    ''' convert a sam to coordinate sorted bam '''
    if (out_bam[-4:] == ".bam"):
        bam_prefix = out_bam[:-4]
    else:
        bam_prefix = out_bam
    bam_unsort = bam_prefix + ".unsorted.bam"
    pysam.view("-Sbu", "-@", str(nthreads), "-o", bam_unsort, in_sam)
    pysam.sort(bam_unsort, "-@", str(nthreads), "-o", bam_prefix + ".bam")
    pysam.index("%s.bam" % bam_prefix)
    os.remove(bam_unsort)


def bam_2_nsort_sam(in_bam, out_sam, nthreads):
    ''' convert a bam to a name sorted sam '''
    bam_sort = in_bam[:-4] + ".nsorted"
    pysam.sort("-n", "-@", str(nthreads),  in_bam, "-o", bam_sort + ".bam")
    pysam.view("-h", "-@", str(nthreads), "-o" , out_sam, bam_sort + ".bam")
    os.remove(bam_sort + ".bam")
