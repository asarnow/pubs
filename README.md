# PUBS 2015 MS Data Analysis

## intensities.py
This script reads MS data and converts it to a convenient form for further analysis.

Other scripts may take advantage using:
`from intensities import *`

Note this also includes `numpy as np`, `scipy as sp` and `pandas as pd`.

Features:

+ Read the proteinGroups.txt table into a PANDAS DataFrames
+ Extract relative intensity matrices for the four sample classes (WCL, WCLp, UB, UBp)
+ Log2-transform experiment : control relative intensity ratios (fold-change)
+ Prepare DataFrames containing the fold-change values and indexed on protein name

## histograms.py
Use intensities.py output to produce histograms of the MS data.

## enrichment.py
Use intensities.py output to compute proteins having sufficiently large
positive or negative fold-change values in a particular experiment.

## global.py
Use itensities.py output for global analysis of all experiments.

Currently:

+ Compute principal components of WCL experiments
+ Plot first 2 PCs' against each other
