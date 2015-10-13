#!/usr/bin/env python2.7
from pylab import *
import pandas
import sys


def main(args):
    # evidence_raw = pandas.read_table('evidence.txt')
    proteins_raw = pandas.read_table('proteinGroups.txt')
    phosphosites_raw = pandas.read_table('Phospho (STY)Sites.txt')

    # evidence = deduplicate(evidence_raw, "Protein Group IDs")
    # evidence = deduplicate(evidence, "Phospho (STY) Site IDs")
    phosphosites = denormalize(phosphosites_raw, "Protein Group IDs")
    evi_norm = pandas.merge(phosphosites,
                            proteins_raw,
                            left_on="Protein Group IDs",
                            right_on="id",
                            how="left")
    evi_norm.to_csv('phospho_norm.txt', sep="\t")


def denormalize(dfin, col):
    df = dfin
    drp = []
    new_entries = []
    for i in reversed(xrange(0, df.shape[0])):
        ids = df.iloc[i][col]
        if type(ids) is not str:
            continue
        tok = ids.split(";")
        if len(tok) > 1:
            idx = [] # needs indices of other columns with multiple values
            # for t in tok:
            for j in xrange(0, len(tok)):
                new_entry = df.iloc[i]
                new_entry[col] = tok[j]
                # for k in idx:
                    # set new_entry columns named in idx to df.iloc[i][k].split(";")[j]
                new_entries.append(new_entry)
            drp.append(i)
    df = df.drop(df.index[drp])
    df = df.append(new_entries)
    return df


if __name__ == "__main__":
    sys.exit(main(sys.argv))
