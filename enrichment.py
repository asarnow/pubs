#!/usr/bin/env python2.7
from intensities import *
import os
import os.path


def threshold_columns(df, col, upper, lower):
    upreg = (df[col] >= upper) & np.isfinite(df[col])
    downreg = (df[col] <= lower) & np.isfinite(df[col])
    return upreg, downreg


def list_genes(enr):
    genes = [t.replace(';', os.linesep) for t in enr]
    genes = [g for g in genes if not (("CON" in g) | ("REV" in g))]
    return genes


def write_enrichment(enr, fname):
    with open(fname, 'w') as f:
        f.write(os.linesep.join(list_genes(enr)) + os.linesep)


def enrichment(df, cols, upper, lower):
    enrichdir = 'data/enrichment'
    for c in df[cols].columns:
        up, down = threshold_columns(df, c, upper, lower)
        write_enrichment(df[up].index, os.path.join(enrichdir, c + '_up' + '.txt'))
        write_enrichment(df[down].index, os.path.join(enrichdir, c + '_down' + '.txt'))


def main():
    enrichment(wcl_foldch, wcl_exp, 3, -3)
    enrichment(wclp_foldch, wclp_exp, 3, -3)
    enrichment(ub_foldch, ub_exp, 3, -3)
    enrichment(ubp_foldch, ubp_exp, 3, -3)

    gsy = pd.read_table('data/go_slim_mapping.tab', header=None, low_memory=False).set_index(0)
    # wcl_go = wcl_foldch.reset_index().merge(
    #         gsy,
    #         how="inner",
    #         left_on="Protein IDs",
    #         right_on=0
    #     ).set_index("Protein IDs")


if __name__ == "__main__":
    main()
