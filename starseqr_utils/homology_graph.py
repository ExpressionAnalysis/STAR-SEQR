#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function)
import os
import logging
import ssw
import numpy as np
import networkx as nx
import pandas as pd
import itertools
import starseqr_utils as su

logger = logging.getLogger('STAR-SEQR')


def get_pairwise_hom(jxn1, jxn2, chim_dir, side):
    hom_res = [0]
    # jxn1
    clean_jxn1 = su.common.safe_jxn(jxn1)
    jxn1_fa = os.path.join(chim_dir, 'transcripts-fusion-' + clean_jxn1 + '.fa')
    # jxn2
    clean_jxn2 = su.common.safe_jxn(jxn2)
    jxn2_fa = os.path.join(chim_dir, 'transcripts-fusion-' + clean_jxn2 + '.fa')
    # get homology scores
    aligner = ssw.Aligner(gap_open=12, gap_extend=4)  # match=2, mismatch=2
    fa1_gen = su.common.fasta_iter(jxn1_fa)
    max_res = 0
    # iterate through transcripts and do SW against fusion partner transcript
    for fa1_head, fa1_seq in fa1_gen:
        fa1_brk = int(fa1_head.split('|')[-1])
        if side == "left":
            fa1_seq = fa1_seq[:fa1_brk]
        else:
            fa1_seq = fa1_seq[fa1_brk:]
        fa2_gen = su.common.fasta_iter(jxn2_fa)
        for fa2_head, fa2_seq in fa2_gen:
            fa2_brk = int(fa2_head.split('|')[-1])
            if side == "left":
                fa2_seq = fa2_seq[:fa2_brk]
            else:
                fa2_seq = fa2_seq[fa2_brk:]
            cmp_align = aligner.align(reference=fa1_seq, query=fa2_seq)
            norm_res = cmp_align.score / (min(len(fa1_seq), len(fa2_seq)) * 2) # normalize hom from 0-1 range
            max_res = max(max_res, norm_res)
            hom_res.append(max_res)
    return(np.max(hom_res))


def prune_homology_graph(df, chim_dir):
    to_remove = []
    df['brk_left_cut'] = df['name'].str.split(":").str[0:3].str.join(sep=":")
    df['brk_right_cut'] = df['name'].str.split(":").str[3:6].str.join(sep=":")
    left_nodes = set(df[df['brk_left_cut'].duplicated()]['brk_left_cut'])
    right_nodes = df[df['brk_right_cut'].duplicated()]['brk_right_cut']
    all_nodes = list(zip(left_nodes, itertools.repeat("left"))) + list(zip(right_nodes, itertools.repeat("right")))
    for node, hom_side in all_nodes:
        node_members = df[((df['brk_' + hom_side + '_cut'] == node))]['name']
        node_graph = nx.Graph()
        node_graph.add_nodes_from(node_members, exprs=10)
        for jxn1, jxn2 in itertools.combinations(node_members, 2):
            pair_score = get_pairwise_hom(jxn1, jxn2, chim_dir, hom_side)
            if pair_score != 0:
                node_graph.add_edge(jxn1, jxn2, weight=pair_score)
        # nx.draw_networkx(node_graph, pos=nx.shell_layout(node_graph), node_size=100)
        # plt.show()
        adj_mat = nx.to_pandas_dataframe(node_graph)
        node_compare = adj_mat[adj_mat.sum()> 0].index.tolist()
        if len(node_compare) > 0:
            node_homdf = df[df['name'].isin(node_compare)][['name', 'TPM_Fusion', 'TPM_Left', 'TPM_Right']].set_index('name')
            node_homdf['max_pairs'] = node_homdf[['TPM_Left','TPM_Right']].max(axis=1)
            node_homdf = node_homdf.sort_values(['TPM_Fusion', 'max_pairs'] , ascending=False)
            node_remove = node_homdf.iloc[1:].index.tolist()
            to_remove.extend(node_remove)
    # use list of to_remove to mark homologous fusions
    return to_remove
