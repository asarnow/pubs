#!/usr/bin/env python2.7
from intensities import *


# 2nd-to-last element is Shmoo / CaCl2.
# Only histogram finite (non-inf, non-NaN) values.
cnts = sp.histogram(
    wcl_foldch["Intensity Shmoo_CaCl2_WCL"][np.isfinite(wcl_foldch["Intensity Shmoo_CaCl2_WCL"])].values)
