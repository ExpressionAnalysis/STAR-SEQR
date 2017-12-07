#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
from six.moves import zip
import six
from six import reraise as raise_
import os
import sys
import re
import gzip
import logging
import time
import errno
import pandas as pd
import numpy as np
from itertools import groupby, repeat
import multiprocessing as mp
import signal
import pysam
from intervaltree_bio import GenomeIntervalTree
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')

try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans


def init_log(name="STAR-SEQR", logfile=None, fileLevel=logging.DEBUG, consoleLevel=logging.INFO):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M:%S")
    trivialformatter = logging.Formatter('%(asctime)-15s - %(levelname)-4s - %(message)s', "%Y-%m-%d %H:%M")

    # create file handler and set level to INFO
    if logfile is None:
        logfile =  sys.exit("No log file defined")
    fh = logging.FileHandler(logfile)
    fh.setLevel(fileLevel)
    fh.setFormatter(trivialformatter)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(consoleLevel)
    ch.setFormatter(trivialformatter)

    # add fh to logger
    logger.addHandler(fh)

    # Show full date time only in header,
    logger.info('\n')
    logger.info('#'*80)
    this_user = str(os.geteuid())
    if not this_user:
        this_user = "USER"
    logger.info('#{0:^78s}#'.format(this_user + time.strftime(" %x  %X")))
    logger.info('#'*80)

    # Now set to normal format
    fh.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def check_file_exists(path):
    if (os.stat(os.path.realpath(path)).st_size == 0):
        logger.error("Exiting. Cannot find file: " + os.path.abspath(path))
        sys.exit(1)


def is_newer_than(target_path, orig_path):
    if not os.path.isfile(target_path):
        return False
    return (os.stat(target_path).st_mtime >= os.stat(orig_path).st_mtime)


def find_resource(filename):
    packagedir = su.__path__[0]
    dirname = os.path.join(packagedir, 'resources')
    fullpath = os.path.abspath(os.path.join(dirname, filename))
    return fullpath


def make_new_dir(newdir):
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            logger.error("Exception: " + str(e))
            raise e
        pass


def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError as e:
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


def pandas_parallel(df, func, nthreads, mp_type, split, *opts):
    '''wrapper to run pandas apply function with multiprocessing. *opts will take any number of optional arguments to be passed
    into the function. Note these will be passed to the function as a tuple and need to be parsed. '''
    def init_worker():
        signal.signal(signal.SIGINT, signal.SIG_IGN)

    logger.info("Begin multiprocessing of function %s in a pool of %s workers using %s protocol" % (str(func.__name__), nthreads, mp_type))
    if split == '':
        logger.debug("*The dataframe will be split evenly across the %s workers" % nthreads)
        df_split = np.array_split(df, min(nthreads, len(df.index)))
    else:
        logger.debug('*The dataframe will be split on the column %s' % split)
        df_split = (df.loc[df[split] == sub, :] for sub in df[split].unique())

    # Pack params for processing
    if not opts:
        params = df_split
    elif len(opts) == 1:
        params = zip(df_split, repeat(opts[0]))
    elif len(opts) > 1:
        params = zip(df_split, repeat(opts))

    init = time.time()
    try:
        logger.debug("*Initializing a %s pool with %s workers" % (mp_type, nthreads))
        pool = mp.Pool(nthreads, init_worker)  # Initialize the pool
        pool_func = getattr(pool, mp_type) # set the function (map, imap, etc)
        if mp_type == "map_async":
            these_res = pool_func(func, params).get()
            out = pd.concat(these_res, axis=0)
        elif mp_type in ["map", "imap", "imap_unordered"]:
            these_res = pool_func(func, params)
            out = pd.concat(these_res, axis=0) # this works even for iterables
        pool.close()
    except KeyboardInterrupt as e:
        logger.error("Error: Keyboard interrupt")
        pool.terminate()
        raise e
    except Exception as e:
        logger.error("Exception: " + str(e))
        pool.terminate()
        traceback = sys.exc_info()[2]
        raise_(ValueError, e, traceback)
    finally:
            pool.join()
    logger.debug("*Time to run pandas_parallel on " + str(func.__name__) + " took %g seconds" % (time.time() - init))
    return(out)


def safe_jxn(jxn):
    ''' make a jxn string safe for folder and filenames'''
    fix_jxn = str(jxn).replace(':', '_')
    fix_jxn = str(fix_jxn).replace('+', 'pos')
    fix_jxn = str(fix_jxn).replace('-', 'neg')
    return fix_jxn


def bed_to_tree(bed):
    with open(bed, 'r') as f:
        btree = GenomeIntervalTree.from_bed(fileobj=f)
    return btree


def subset_bed_func(jxn, targets_tree, sub_style='either'):
    chrom1, pos1, str1, chrom2, pos2, str2, repleft, repright = re.split(':', jxn)
    intersect1 = len(targets_tree[str(chrom1)].search(int(pos1)))
    intersect2 = len(targets_tree[str(chrom2)].search(int(pos2)))
    if sub_style == 'either':
        if intersect1 or intersect2 >= 1:
            return 1
        else:
            return 0
    elif sub_style == 'both':
        if intersect1 >= 1 and intersect2 >= 1:
            return 1
        else:
            return 0


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def fasta_iter(fasta_name):
    ''' Given a fasta file, yield tuples of header, sequence'''
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
    for header in faiter:
        header = six.next(header)[1:].strip()  # drop the '>'
        # join all sequence lines to one.
        seq = ''.join(s.strip() for s in six.next(faiter))
        yield header, seq
    fh.close()


class FastqRead(object):
    """Represents 1 read, or 4 lines of a casava 1.8-style fastq file."""
    def __init__(self, file_obj):
        header = six.next(file_obj).rstrip()
        assert header.startswith('@')
        self.header = header
        self.sequence = six.next(file_obj).rstrip()
        self.seq_len = len(self.sequence)
        six.next(file_obj)
        self.quality = six.next(file_obj).rstrip()

    def __str__(self):
        return "{0}\n{1}\n+\n{2}\n+\n{3}\n".format(self.header, self.sequence,
                                                   self.quality, self.seq_len)


class FastqParser(object):
    """Parses a fastq file into FastqReads."""
    def __init__(self, filename):
        if filename.endswith('.gz'):
            self._file = gzip.open(filename, 'rb')
        else:
            self._file = open(filename, 'r')
        self._line_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        read = FastqRead(self._file)
        if read:
            self._line_index += 4
            return read
        else:
            self._file.close()
            raise StopIteration()

    next = __next__  # Python 2


def sam_2_coord_bam(in_sam, out_bam, nthreads):
    ''' convert a sam to coordinate sorted bam '''
    if (out_bam[-4:] == ".bam"):
        bam_prefix = str(out_bam[:-4])
    else:
        bam_prefix = str(out_bam)
    bam_unsort = bam_prefix + ".unsorted.bam"
    pysam.view("-Sbu", "-@", str(nthreads), "-o", str(bam_unsort), str(in_sam), catch_stdout=False)
    pysam.sort(str(bam_unsort), "-@", str(nthreads), "-o", str(bam_prefix + ".bam"), catch_stdout=False)
    pysam.index("%s.bam" % bam_prefix, catch_stdout=False)
    os.remove(bam_unsort)


def bam_2_nsort_sam(in_bam, out_sam, nthreads):
    ''' convert a bam to a name sorted sam '''
    bam_sort = str(in_bam[:-4] + ".nsorted")
    pysam.sort("-n", "-@", str(nthreads), str(in_bam), "-o", str(bam_sort + ".bam"), catch_stdout=False)
    pysam.view("-h", "-@", str(nthreads), "-o", str(out_sam), str(bam_sort + ".bam"), catch_stdout=False)
    os.remove(bam_sort + ".bam")


def sortnindex_bam(in_bam, out_bam, nthreads):
    pysam.sort("-@", str(nthreads), str(in_bam), "-o", str(out_bam), catch_stdout=False)
    pysam.index(out_bam, catch_stdout=False)
    os.remove(in_bam)


def index_bam(in_bam):
    if not os.path.isfile(in_bam + ".bai") or os.stat(in_bam).st_mtime >= os.stat(in_bam + ".bai").st_mtime:
        pysam.index(in_bam, catch_stdout=False)
