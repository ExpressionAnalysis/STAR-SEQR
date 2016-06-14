#!/usr/bin/env python

from __future__ import print_function
import sys
import time
import math
import pandas as pd
import numpy as np
import logging


logger = logging.getLogger("STAR-SEQR")


def write_header_vcf(args, fh):
    '''Follow the spec here: https://samtools.github.io/hts-specs/VCFv4.2.pdf'''
    today = time.strftime('%m/%d/%Y')
    vpcmd = str(' '.join(sys.argv))
    data_head = '\t'.join(
        ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.prefix])
    header = '\n'.join(['##fileformat=VCFv4.2', '##fileDate=' + today,
                        '##source=STAR-SV', '##reference=' + 'hg19',
                        '##FILTER=<ID=PASS,Description="All filters passed">',
                        '##FILTER=<ID=FAIL,Description="Site failed to reach confidence">',
                        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
                        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
                        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',
                        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
                        '##INFO=<ID=HOMLEN,Number=-1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">',
                        '##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">',
                        '##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
                        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
                        '##ALT=<ID=DEL,Description="Deletion">',
                        '##ALT=<ID=DUP,Description="Duplication">',
                        '##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">',
                        '##ALT=<ID=INS,Description="Insertion of novel sequence">',
                        '##ALT=<ID=INV,Description="Inversion">',
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                        '##FORMAT=<ID=GQ,Number=1,Type=String,Description="Genotype Quality">',
                        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
                        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth at position">',
                        '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depths of ref and alt alleles in that order">',
                        '##FORMAT=<ID=FSAD,Number=2,Type=Integer,Description="Allele depths on forward strands of ref and alt alleles in that order">',
                        '##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias">',
                        '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average Base Quality of ALT">',
                        '##FORMAT=<ID=PV,Number=1,Type=Float,Description="Adjusted p-value">',
                        '##variants_justified=left',
                        '##STAR-SV_CMD=' + vpcmd,
                        data_head
                        ])
    print(header, file=fh)


def process(svtable, args):
    vcf_sv_fh = open(args.prefix + '_sv.vcf', 'w')
    write_header_vcf(args, vcf_sv_fh)
