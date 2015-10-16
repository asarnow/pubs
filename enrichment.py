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


def write_enrichment(enr, fname):
    with open(os.path.join(enrichdir, fname), 'w') as f:
        f.write('\n'.join(list_genes(enr)))


wcl_cacl2_up, wcl_cacl2_down = enrich(wcl_foldch, "Intensity Shmoo_CaCl2_WCL", 3, -3)
wcl_cmk1_up, wcl_cmk1_down = enrich(wcl_foldch, "Intensity Shmoo_Cmk1KO_WCL", 3, -3)
wclp_cmk1_up, wclp_cmk1_down = enrich(wclp_foldch, "Intensity Shmoo_Cmk1KO_WCLP", 3, -3)

list_genes = lambda x: [t.replace(';', '\n') for t in x]

write_enrichment(wcl_cacl2_up, 'WCL_CaCl2_up_gt3.txt')
write_enrichment(wcl_cacl2_down, 'WCL_CaCl2_down_le3.txt')
write_enrichment(wcl_cmk1_up, 'WCL_CaCl2_up_gt3.txt')
write_enrichment(wcl_cmk1_down, 'WCL_CaCl2_down_le3.txt')
write_enrichment(wclp_cmk1_up, 'WCLp_CaCl2_up_gt3.txt')
write_enrichment(wclp_cmk1_down, 'WCLp_CaCl2_down_le3.txt')
