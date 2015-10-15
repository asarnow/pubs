#!/usr/bin/env python2.7
from intensities import *


def enrich(df, col, upper, lower):
    upreg = df[
        (df[col] > upper) &
        np.isfinite(df[col])].index
    downreg = df[
        (df[col] < lower) &
        np.isfinite(df[col])].index
    return upreg, downreg

# target = wcl_foldch[
# ((wcl_foldch["Intensity Shmoo_CaCl2_WCL"] > 3) |
# (wcl_foldch["Intensity Shmoo_CaCl2_WCL"] < -3)) &
# np.isfinite(wcl_foldch["Intensity Shmoo_CaCl2_WCL"])].index
# target = wcl_foldch[
#     (wcl_foldch["Intensity Shmoo_CaCl2_WCL"] > 3) &
#     np.isfinite(wcl_foldch["Intensity Shmoo_CaCl2_WCL"])].index

wcl_cacl2_up, wcl_cacl2_down = enrich(wcl_foldch, "Intensity Shmoo_CaCl2_WCL", 3, -3)

list_genes = lambda x: [t.replace(';', '\n') for t in x]

with open('data/enrichment/WCL_CaCl2_up_gt3.txt', 'w') as f:
    f.write('\n'.join(list_genes(wcl_cacl2_up)))
