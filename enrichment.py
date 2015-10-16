#!/usr/bin/env python2.7
from intensities import *
import os.path


def enrich(df, col, upper, lower):
    upreg = df[
        (df[col] >= upper) &
        np.isfinite(df[col])].index
    downreg = df[
        (df[col] <= lower) &
        np.isfinite(df[col])].index
    return upreg, downreg


enrichdir = 'data/enrichment'


def list_genes(enr):
    genes = [t.replace(';', '\n') for t in enr]
    genes = [g for g in genes if not (("CON" in g) | ("REV" in g))]
    return genes


def write_enrichment(enr, fname):
    with open(fname, 'w') as f:
        f.write('\n'.join(list_genes(enr)) + '\n')


def enrichment(df, cols, upper, lower):
    for c in df[cols].columns:
        up, down = enrich(df, c, upper, lower)
        write_enrichment(up, os.path.join(enrichdir, c + '_up' + '.txt'))
        write_enrichment(down, os.path.join(enrichdir, c + '_down' + '.txt'))


enrichment(wcl_foldch, wcl_exp, 3, -3)
enrichment(wclp_foldch, wclp_exp, 3, -3)
enrichment(ub_foldch, ub_exp, 3, -3)
enrichment(ubp_foldch, ubp_exp, 3, -3)

gsy = pd.read_table('data/go_slim_mapping.tab', low_memory=False)
