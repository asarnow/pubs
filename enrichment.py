#!/usr/bin/env python2.7
from intensities import *
import sys
import os
import os.path
import go_obo_parser as obo
import scipy.stats as stats


def threshold_columns(df, col, upper=3., lower=-3.):
    upreg = (df[col] >= upper) & np.isfinite(df[col])
    downreg = (df[col] <= lower) & np.isfinite(df[col])
    return upreg, downreg


def list_genes(enr):
    genes = [t.replace(';', os.linesep) for t in enr]
    genes = [g for g in genes if not (("CON" in g) | ("REV" in g))]
    return genes


def write_genes(enr, fname):
    with open(fname, 'w') as f:
        f.write(os.linesep.join(list_genes(enr)) + os.linesep)


def threshold_and_write(df, cols, upper, lower):
    enrichdir = 'data/threshold'
    for c in df[cols].columns:
        up, down = threshold_columns(df, c, upper, lower)
        write_genes(df[up].index, os.path.join(enrichdir, c + '_up' + '.txt'))
        write_genes(df[down].index, os.path.join(enrichdir, c + '_down' + '.txt'))


def enrich(df, upper=3., lower=-3., bonferroni=False):
    gsy = pd.read_table('data/go_slim_mapping.tab', header=None, low_memory=False).set_index(0)
    goslim = [g for g in obo.parse_go_obo('data/goslim_yeast.obo')]
    if bonferroni:
        bcor = len(set(gsy[5]))
    else:
        bcor = 1.
    # wcl_go = wcl_foldch.reset_index().merge(
    #     gsy, how="inner", left_on="Protein IDs", right_on=0).set_index("Protein IDs")
    res_up = {c: [] for c in df.columns}
    res_down = {c: [] for c in df.columns}
    for c in df.columns:
        up, down = threshold_columns(df, c, upper, lower)
        goids_up = gsy.loc[df.index[up]][5]
        goids_down = gsy.loc[df.index[down]][5]
        goids_global = gsy.loc[df.index][5].tolist()
        terms_up = [g for g in goslim if g['id'] in set(goids_up)]
        terms_down = [g for g in goslim if g['id'] in set(goids_up)]
        for tpl in calc_enrichment(terms_up, goids_up, goids_global, bcor):
            res_up[c].append(tpl)
        for tpl in calc_enrichment(terms_down, goids_down, goids_global, bcor):
            res_down[c].append(tpl)
    return res_up, res_down


def calc_enrichment(terms, goids, goids_global, bcor):
    goidsr = goids.reset_index()
    goidsl = goidsr[5].tolist()
    for t in terms:
        cnt = goidsl.count(t['id'])
        global_cnt = goids_global.count(t['id'])
        pval = stats.hypergeom.sf(cnt, len(goids_global), global_cnt, len(goids)) * bcor
        if pval < 0.05:
            yield (t['id'], t['name'], pval, ', '.join(goidsr[0][goidsr[5] == t['id']]))


def write_enrichment(enr, fname):
    with open(fname, 'w') as f:
        for tpl in enr:
            f.write('\t'.join(str(s) for s in tpl) + os.linesep)


def calc_and_write(df, upper=3., lower=-3., bonferroni=False):
    enr_up, enr_down = enrich(df, upper, lower, bonferroni)
    enrichdir = 'data/enrichment'
    for k in enr_up:
        fn = os.path.join(enrichdir, k + "_up_go.txt")
        write_enrichment(enr_up[k], fn)
    for k in enr_down:
        fn = os.path.join(enrichdir, k + "_down_go.txt")
        write_enrichment(enr_down[k], fn)


def main(args):
    if args.thresh:
        threshold_and_write(wcl_foldch, wcl_exp, args.upper, args.lower)
        threshold_and_write(wclp_foldch, wclp_exp, args.upper, args.lower)
        threshold_and_write(ub_foldch, ub_exp, args.upper, args.lower)
        threshold_and_write(ubp_foldch, ubp_exp, args.upper, args.lower)
    if args.go:
        calc_and_write(wcl_foldch, args.upper, args.lower, args.bonferroni)
        calc_and_write(wclp_foldch, args.upper, args.lower, args.bonferroni)
        calc_and_write(ub_foldch, args.upper, args.lower, args.bonferroni)
        calc_and_write(ubp_foldch, args.upper, args.lower, args.bonferroni)
    else:
        return 1
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Threshold protein abundances and compute GO enrichment scores.")
    parser.add_argument("-u", "--upper", action="store", type=float, default=3., dest="upper")
    parser.add_argument("-l", "--lower", action="store", type=float, default=-3., dest="lower")
    parser.add_argument("-t", "--thresholds", action="store_true", default=False, dest="thresh")
    parser.add_argument("-g", "--go", action="store_true", default=False, dest="go")
    parser.add_argument("-b", "--bonferroni", action="store_true", default=False, dest="bonferroni")
    sys.exit(main(parser.parse_args()))
