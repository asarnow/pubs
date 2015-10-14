# PUBS 2015 MS Data Analysis

## msanal.py
msanal.py is a skeletal MS data analysis program.

Currently, it supports the following operations:

+ Parse data files using PANDAS
+ "Denormalize" rows containing fields with multiple values
+ Use table join to merge multiple MS data files based on a shared field

## intensities.py
intensities.py is another skeletal MS data analysis script.
It's currently written to run in ipython after `%pylab` magic.
Suggest copying the whole code and using `%paste` magic.

+ Read the proteinGroups.txt table into a PANDAS data frame
+ Extract relative intensity matrices for the four sample classes (WCL, WCLp, UB, UBp)
+ Log2-transform experiment : control relative intensity ratios (fold-change)