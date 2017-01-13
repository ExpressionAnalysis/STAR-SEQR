#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function)
import sys
from argparse import ArgumentParser


def create_bedpe(file_in, file_out):
    with open(file_out, 'w') as file_out_fh:
        with open(file_in, 'r') as starOutput:
            for line in starOutput:
                if not line.startswith(('#','NAME')):
                    vals = line.strip().split()
                    valsL = vals[6].split(':')
                    valsR = vals[7].split(':')
                    chrL = valsL[0]
                    posL = valsL[1]
                    geneL = vals[0].split('--')[0]
                    strandL = valsL[2]
                    chrR = valsR[0]
                    posR = valsR[1]
                    geneR = vals[0].split('--')[1]
                    strandR = valsR[2]
                    total_uniqreads = int(vals[1]) + int(vals[2]) + int(vals[3])
                    bedpe = '\t'.join([chrL, str(int(posL)-1), posL, chrR, str(int(posR)-1),
                                       posR,'-'.join([geneL, geneR]), str(total_uniqreads), strandL, strandR])
                    print(bedpe, file=file_out_fh)


if __name__ == '__main__':
    # get args from command line
    usage = 'starseqr2bedpe.py -i file.breakpoints.txt -o file.bedpe'
    parser = ArgumentParser(description="Convert breakpoints.txt to a bedpe", epilog=usage)
    parser.add_argument('-i', '--input', type=str, required=True, help='breakpoints.txt input')
    parser.add_argument('-o', '--output', type=str, required=True, help='bedpe out file, /dev/stdout will work.')
    args = parser.parse_args()
    # run the function
    create_bedpe(args.input, args.output)