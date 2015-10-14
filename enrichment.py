#!/usr/bin/env python2.7
from intensities import *

# target = wcl_foldch[((wcl_foldch["Intensity Shmoo_CaCl2_WCL"] > 3) | (wcl_foldch["Intensity Shmoo_CaCl2_WCL"] < -3)) & np.isfinite(wcl_foldch["Intensity Shmoo_CaCl2_WCL"])].index
target = wcl_foldch[
    (wcl_foldch["Intensity Shmoo_CaCl2_WCL"] > 3) & np.isfinite(wcl_foldch["Intensity Shmoo_CaCl2_WCL"])].index

target = [t.replace(';', '\n') for t in target]
with open('data/enrichment/WCL_CaCl2_up_gt3.txt', 'w') as f:
    f.write('\n'.join(target))
