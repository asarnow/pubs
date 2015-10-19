#!/usr/bin/env python2.7
from intensities import *
import matplotlib.pyplot as plt
import seaborn


# 2nd-to-last element is Shmoo / CaCl2.
# Only histogram finite (non-inf, non-NaN) values.

plt.hist(ub_foldch['Intensity Shmoo_Cmk1KO_Ub'][np.isfinite(ub_foldch['Intensity Shmoo_Cmk1KO_Ub'])].values,50)
plt.axvline(x=2.0, ymin= 0, linewidth=2, color='r')
plt.axvline(x=-2.0, ymin= 0, linewidth=2, color='r')
plt.ylabel('Number of Proteins')
plt.xlabel('log2 Fold Change to WT Experiment')
plt.title('CMK1 KO Ub FC to WT')
plt.grid(False)
plt.show()
