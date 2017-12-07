#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import sys
import logging
import subprocess as sp
import starseqr_utils as su

logger = logging.getLogger("STAR-SEQR")


def run_star(prefix, fq1, fq2, star_index, threads, mode):
    '''
    Run STAR alignment for RNA or DNA with different sensitivity parameters.
    '''
    logger.info("Starting STAR Alignment")
    if not os.path.isfile(prefix + ".Chimeric.out.junction"):
        if not os.path.isdir(star_index):
            logger.error("Error: STAR index was not found at " + star_index + " Please update the path.", exc_info=True)
            sys.exit(1)
        # outFilterMultimapNmax does not have any effect on chimeric alignments
        # chimSegmentReadGapMax allows for variation at breakpoints
        # chimScoreDropMax is the difference in aligned read length and total read length due to clipping
        # chimMainSegmentMultNmax - new parameter >2.5b to remove chimeric multimappers
        # withinBam and SeparateSAMold give equivalent results
        # mismatches at sj: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.
        STAR_args = ['STAR', '--readFilesIn', fq1, fq2, '--readFilesCommand', 'zcat',
                     '--runThreadN', str(threads), '--genomeDir', star_index,
                     '--outFileNamePrefix ', prefix + ".", '--chimScoreJunctionNonGTAG', -1,
                     '--outSAMtype', 'None', '--chimOutType', 'SeparateSAMold',
                     # '--outSAMtype', 'BAM', 'SortedByCoordinate', '--chimOutType', 'WithinBAM',
                     '--alignSJDBoverhangMin', 5, '--outFilterMultimapScoreRange', 1,
                     '--outFilterMultimapNmax', 5,
                     '--outMultimapperOrder', 'Random', '--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'ch']
        # choose sensitivity mode
        if (mode == 0):
            sens_params = ['--chimSegmentMin', 10, '--chimJunctionOverhangMin', 10,
                           '--chimScoreMin', 1, '--chimScoreDropMax', 20,
                           '--chimScoreSeparation', 10, '--chimSegmentReadGapMax', 3,
                           '--chimFilter', 'None', '--twopassMode', "None",
                           '--alignSJstitchMismatchNmax', 5, -1, 5, 5,
                           '--chimMainSegmentMultNmax', 1]
        elif (mode == 1):
            sens_params = ['--chimSegmentMin', 10, '--chimJunctionOverhangMin', 10,
                           '--chimScoreMin', 1, '--chimScoreDropMax', 30,
                           '--chimScoreSeparation', 7, '--chimSegmentReadGapMax', 3,
                           '--chimFilter', 'None', '--twopassMode', "None",
                           '--alignSJstitchMismatchNmax', 5, -1, 5, 5,
                           '--chimMainSegmentMultNmax', 10]
        STAR_args.extend(sens_params)
        # Need to convert all to string
        STAR_args = list(map(str, STAR_args))

        logger.info("*STAR Command: " + " ".join(STAR_args))
        # run STAR
        try:
            p = sp.Popen(STAR_args, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()
            if stdout:
                logger.info(stdout)
            if stderr:
                logger.error(stderr)
            if p.returncode != 0:
                logger.error("Error: STAR failed", exc_info=True)
                sys.exit(1)
        except (OSError) as o:
            logger.error("Exception: " + str(o))
            logger.error("STAR Failed", exc_info=True)
            sys.exit(1)
        su.common.check_file_exists(prefix + ".Chimeric.out.junction")
        su.common.check_file_exists(prefix + ".Chimeric.out.sam")
    else:
        su.common.check_file_exists(prefix + ".Chimeric.out.junction")
        su.common.check_file_exists(prefix + ".Chimeric.out.sam")
        logger.warn("Skipping STAR alignment as files already exist!")
    logger.info("STAR Alignment Finished!")
