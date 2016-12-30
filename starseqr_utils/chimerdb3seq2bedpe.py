#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function)
import sys
import time
from argparse import ArgumentParser
import pandas as pd


def parse_args():
    usage = " "
    parser = ArgumentParser(
        description="convert FusionCancer DB to bedpe:", epilog=usage)
    parser.add_argument('-i', '--file', type=str, required=True,
                        help='excel file from chimerdb3')
    parser.add_argument('-i2', '--file2', type=str, required=True,
                        help='excel file from chimerdb3 with fixed coordinates for chimerdb2 db')
    # parser.add_argument('-k', '--file3', type=str, required=True,
    #                     help='excel file from chimerdb3 for KB database')
    parser.add_argument('-o', '--out', type=str, required=False,
                        help='name of out file')
    return parser.parse_args()


def main():
    start = time.time()
    args = parse_args()
    chimerdb_seq = pd.read_excel(args.file, header=0)
    chimerdb_seq = chimerdb_seq[chimerdb_seq['Genome_Build_Version'] == 'hg19']
    # update names
    chimerdb_seq['Source'] = chimerdb_seq['Source'].str.replace('PRADA', 'TCGA-PRADA')
    chimerdb_seq['Source'] = chimerdb_seq['Source'].str.replace('FusionScan', 'TCGA-FusionScan')
    chimerdb_seq['Source'] = chimerdb_seq['Source'].str.replace('TopHat-Fusion', 'TCGA-TopHat-Fusion')
    # file2
    chimerdb_seq2 = pd.read_excel(args.file2, header=0)
    chimerdb_seq2 = chimerdb_seq2[chimerdb_seq2['Genome_Build_Version'] == 'hg19']
    chimerdb_all = pd.concat([chimerdb_seq, chimerdb_seq2])
    chimerdb_all['Cancer'] = chimerdb_all['Cancertype_or_disease'].str.replace(' ','_')
    # aggregate by left chr:pos and right chr:pos
    chimerdb_all.update(chimerdb_all[['ChimerKB', 'ChimerPub','Kinase' , 'Oncogene','Tumor_suppressor',
                                      'Receptor', 'Transcription_Factor', 'Cancer']].fillna('NA'))
    grouped_chimerdb = chimerdb_all.groupby(['H_chr', 'H_position', 'T_chr', 'T_position', 'Fusion_pair',
                                          'ChimerKB', 'ChimerPub','Kinase' , 'Oncogene','Tumor_suppressor',
                                          'Receptor', 'Transcription_Factor',
                                         ], as_index=True)
    df = grouped_chimerdb.agg({'Source': {'Source': lambda col: ','.join(set(col))},
                               'Cancer': {'Cancer': lambda col: ','.join(set(col))},
                               'BarcodeID': 'count'}).reset_index()
    df.columns = ['H_chr', 'H_position', 'T_chr', 'T_position', 'Fusion_pair', 'ChimerKB', 'ChimerPub',
                  'Kinase', 'Oncogene', 'Tumor_suppressor', 'Receptor', 'Transcription_Factor',
                  'Source', 'Counts', 'Cancer']
    # merge features into a single string
    df['Features'] = df[['Kinase', 'Oncogene', 'Tumor_suppressor', 'Receptor',
                         'Transcription_Factor', 'ChimerKB',
                         'ChimerPub']].apply(lambda x: ','.join([str(y) for y in x if y != 'NA']), axis=1)
    df.replace('', 'NA', inplace=True)

    # get columns for bedpe format
    df['c1'] = df['H_chr']
    df['c2'] = df['T_chr']
    df['s1'] = df['H_position']
    df['s2'] = df['T_position']
    df['s1'] = df['s1'].astype('int')
    df['s2'] = df['s2'].astype('int')
    df['e1'] = df['s1'] + 1
    df['e2'] = df['s2'] + 1
    df['strand1'] = '.'
    df['strand2'] = '.'
    df.sort_values(['c1', 's1'], ascending=[True, True], inplace=True)
    outcols = ['c1', 's1', 'e1', 'c2', 's2', 'e2', 'Fusion_pair', 'Counts', 'strand1', 'strand2', 'Source', 'Cancer', 'Features']

    df.to_csv(path_or_buf=args.out, header=False, sep='\t', columns=outcols, mode='w', index=False)

    end = time.time()
    elapsed = end - start
    print("Program took  %g seconds" % (elapsed))


if __name__ == "__main__":
    sys.exit(main())