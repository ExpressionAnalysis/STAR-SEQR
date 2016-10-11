#!/usr/bin/env python

from __future__ import print_function
import sys
import time
from argparse import ArgumentParser
import pandas as pd


def parse_args():
    usage = " "
    parser = ArgumentParser(
        description="convert FusionCancer DB to bedpe:", epilog=usage)
    parser.add_argument('-i', '--file', type=str, required=True,
                        help='excel file from FusionCancer')
    parser.add_argument('-o', '--out', type=str, required=False,
                        help='name of out file')
    return parser.parse_args()


def main():
    start = time.time()
    args = parse_args()
    df = pd.read_excel(args.file, header=0)

    df['c1'],  df['s1'] = zip(*df['Head_breakpoint_location'].str.split(':').tolist())
    df['c2'],  df['s2'] = zip(*df['Tail_breakpoint_location'].str.split(':').tolist())
    df['s1'] = df['s1'].astype('int')
    df['s2'] = df['s2'].astype('int')
    df['e1'] = df['s1'] + 1
    df['e2'] = df['s2'] + 1
    outcols = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Database_ID', 'Rate', 'NA', 'NA', 'Fusion_type', 'Tool_n']
    df.to_csv(path_or_buf=args.out, header=False, sep='\t', columns=outcols, mode='w', index=False)

    end = time.time()
    elapsed = end - start
    print("Program took  %g seconds" % (elapsed))


if __name__ == "__main__":
    sys.exit(main())