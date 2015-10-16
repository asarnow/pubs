#!/usr/bin/env python2.7
from intensities import *
import os.path


def enrich(df, col, upper, lower):
    upreg = df[
        (df[col] > upper) &
        np.isfinite(df[col])].index
    downreg = df[
        (df[col] < lower) &
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


for c in wcl_foldch[wcl_exp].columns:
    up, down = enrich(wcl_foldch, c, 3, -3)
    write_enrichment(up, os.path.join(enrichdir, c + '.txt'))
