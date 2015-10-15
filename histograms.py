#!/usr/bin/env python2.7
from intensities import *


# 2nd-to-last element is Shmoo / CaCl2.
# Only histogram finite (non-inf, non-NaN) values.
cnts = sp.histogram(wcl_foldch[wcl_foldch.columns[-2]][np.isfinite(wcl_foldch[wcl_foldch.columns[-2]])].values)
