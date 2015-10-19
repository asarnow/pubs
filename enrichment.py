#!/usr/bin/env python2.7
from intensities import *
import sys
import os
import os.path
import copy
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


def threshold_and_write(enrichdir, df, cols, upper, lower):
    for c in df[cols].columns:
        up, down = threshold_columns(df, c, upper, lower)
        write_genes(df[up].index, os.path.join(enrichdir, c + '_up' + '.txt'))
        write_genes(df[down].index, os.path.join(enrichdir, c + '_down' + '.txt'))


def enrich(df, upper=3., lower=-3., alpha=0.05, correction=None):
    gsy = pd.read_table('data/go_slim_mapping.tab', header=None, low_memory=False).set_index(0)
    goslim = [g for g in obo.parse_go_obo('data/goslim_yeast.obo')]
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
        for lst in calc_enrichment(terms_up, goids_up, goids_global):
            res_up[c].append(lst)
        for lst in calc_enrichment(terms_down, goids_down, goids_global):
            res_down[c].append(lst)
        if correction == "holm":
            res_up[c] = holm(res_up[c], alpha)
            res_down[c] = holm(res_down[c], alpha)
        elif correction == "bonferroni":
            res_up[c] = bonferroni(res_up[c], len(set(gsy[5])), alpha)
            res_down[c] = bonferroni(res_down[c], len(set(gsy[5])), alpha)
        else:
            res_up[c] = [lst for lst in res_up[c] if lst[2] < alpha]
            res_down[c] = [lst for lst in res_down[c] if lst[2] < alpha]
    return res_up, res_down


def calc_enrichment(terms, goids, goids_global):
    goidsr = goids.reset_index()
    goidsl = goidsr[5].tolist()
    for t in terms:
        cnt = goidsl.count(t['id'])
        global_cnt = goids_global.count(t['id'])
        pval = stats.hypergeom.sf(cnt, len(goids_global), global_cnt, len(goids))
        yield [t['id'], t['name'], pval, ', '.join(goidsr[0][goidsr[5] == t['id']])]


def write_enrichment(enr, fname):
    with open(fname, 'w') as f:
        for tpl in enr:
            f.write('\t'.join(str(s) for s in tpl) + os.linesep)


def calc_and_write(enrichdir, df, upper=3., lower=-3., alpha=0.05, correction=None):
    enr_up, enr_down = enrich(df, upper, lower, alpha, correction)
    for k in enr_up:
        fn = os.path.join(enrichdir, k + "_up_go.txt")
        write_enrichment(enr_up[k], fn)
    for k in enr_down:
        fn = os.path.join(enrichdir, k + "_down_go.txt")
        write_enrichment(enr_down[k], fn)


def bonferroni(res, bcor, alpha):
    new = []
    for lst in res:
        if lst[2] * bcor < alpha:
            kst = copy.deepcopy(lst)
            kst[2] = lst[2] * bcor
            new.append(kst)
    return new


def holm(res, alpha):
    pvs = sorted(res, key=lambda x: x[2])
    k = 0
    for i in xrange(0, len(pvs)):
        if pvs[i][2] > alpha / (len(pvs) + 1 - i):
            k = i
            break
    return pvs[0:k]


def main(args):
    if args.thresh:
        threshold_and_write(args.dir, wcl_foldch, wcl_exp, args.upper, args.lower)
        threshold_and_write(args.dir, wclp_foldch, wclp_exp, args.upper, args.lower)
        threshold_and_write(args.dir, ub_foldch, ub_exp, args.upper, args.lower)
        threshold_and_write(args.dir, ubp_foldch, ubp_exp, args.upper, args.lower)
    if args.go:
        if args.bonferroni:
            cor = "bonferroni"
        elif args.holm:
            cor = "holm"
        else:
            cor = None
        calc_and_write(args.dir, wcl_foldch, args.upper, args.lower, args.alpha, cor)
        calc_and_write(args.dir, wclp_foldch, args.upper, args.lower, args.alpha, cor)
        calc_and_write(args.dir, ub_foldch, args.upper, args.lower, args.alpha, cor)
        calc_and_write(args.dir, ubp_foldch, args.upper, args.lower, args.alpha, cor)
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
    parser.add_argument("-s", "--holm", action="store_true", default=False, dest="holm")
    parser.add_argument("-a", "--alpha", action="store", type=float, default=0.05, dest="alpha")
    parser.add_argument("-d", "--dir", action="store", default=".", dest="dir")
    sys.exit(main(parser.parse_args()))
