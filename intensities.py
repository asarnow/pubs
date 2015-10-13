import pandas


proteins = pandas.read_table('data/pubs2015/proteinGroups.txt', low_memory=False)
intensity_cols = [c for c in proteins.columns if 'intensity ' in c.lower() and 'lfq' not in c.lower()]
ub_cols = [c for c in intensity_cols if '_ub' in c.lower()]
wcl_cols = [c for c in intensity_cols if '_wcl' in c.lower()]
ubp_cols = [c for c in intensity_cols if '_ubp' in c.lower()]
wclp_cols = [c for c in intensity_cols if '_wclp' in c.lower()]
mask = (proteins['Reverse'] != '+') & (proteins['Potential contaminant'] != '+')
intensities = proteins[mask][intensity_cols]
total_intensities = proteins[intensity_cols].sum(axis=0)
normed_intensities = intensities / total_intensities
idx = (normed_intensities != 0).any(axis=1)
names = proteins[mask][idx]['Protein IDs']
nonzero_intensities = normed_intensities[idx]

wcl = nonzero_intensities[wcl_cols]
wclp = nonzero_intensities[wclp_cols]
ub = nonzero_intensities[ub_cols]
ubp = nonzero_intensities[ubp_cols]


