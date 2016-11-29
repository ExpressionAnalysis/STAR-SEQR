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

    df['Fusion_name'] = df['Head_gene_symbol'].astype(str) + "--" + df['Tail_gene_symbol'].astype(str)
    df['Cancer_fix'] = df['Cancer_type'].str.replace(' ','_')
    grouped_fuca = df.groupby([ 'c1', 's1', 'c2', 's2', 'Database_ID',], as_index=True)
    fuca = grouped_fuca.agg({'Cancer_fix': {'Cancer_fix': lambda col: ','.join(set(col))},
                               'Detected_tools': {'Detected_tools': lambda col: ','.join(set(col))},
                               'Fusion_name': {'Fusion_name': lambda col: ','.join(set(col))},
                               'Rate': 'max'}).reset_index()
    fuca.columns = ['c1', 's1', 'c2', 's2', 'Database_ID', 'Detected_tools', 'Rate', 'Cancer_fix', 'Fusion_name']
    fuca['Methods'] = fuca['Detected_tools'].str.split(",", expand=False).apply(lambda x: ','.join(set(x)))
    fuca['Cancers'] = fuca['Cancer_fix'].str.split(",", expand=False).apply(lambda x: ','.join(set(x)))
    fuca.sort_values(['c1', 's1'], ascending=[True, True], inplace=True)

    fuca['s1'] = fuca['s1'].astype('int')
    fuca['s2'] = fuca['s2'].astype('int')
    fuca['e1'] = fuca['s1'] + 1
    fuca['e2'] = fuca['s2'] + 1
    fuca['strand1'] = '.'
    fuca['strand2'] = '.'

    outcols = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Database_ID', 'Rate', 'strand1', 'strand2', 'Methods', 'Cancers', 'Fusion_name']
    fuca.to_csv(path_or_buf=args.out, header=False, sep='\t', columns=outcols, mode='w', index=False)

    end = time.time()
    elapsed = end - start
    print("Program took  %g seconds" % (elapsed))


if __name__ == "__main__":
    sys.exit(main())