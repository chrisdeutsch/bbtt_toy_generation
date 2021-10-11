#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

from utils import masspoints, addHeavyFlavourSplit
from hadhad_utils import getHadhadDf, edgesHadhad


parser = argparse.ArgumentParser()
parser.add_argument("ntuple")
parser.add_argument("-o", "--outfile", required=True)
args = parser.parse_args()


treenames = [
    "data",
    "Htautau", "VH",
    "ttH", "ttV",
    "Wenu", "Wmunu", "Wtaunu",
    "Zee", "Zmumu", "Ztautau",
    "Diboson",
    "Fake",
    "singletop", "ttbar", "ttbarFakesMC",
]

dfs = []
for tree in treenames:
    dfs.append(getHadhadDf(args.ntuple, tree))

df = pd.concat(dfs)
del dfs

# Adding split of Z+jets in bb, bc, cc, bl, cl, l
addHeavyFlavourSplit(df)
df["sample"] = df["sample"].astype("category")

# Remove Zttl (also removed in workspaces)
df = df.loc[df["sample"] != "Zttl"]

# Yield tables
print("Yield table")
print(
    df.groupby("sample")["weight"].agg(
        Entries="count",
        Integral="sum")
)

# Adding bin indices
for mass in masspoints:
    # Get bin edges
    edges = edgesHadhad[mass]

    # Discretize MVA scores according to binning
    idx = np.digitize(df[f"PNN{mass}"], bins=edges)

    # Fits uint8
    assert np.all((idx >= 0) & (idx <= 255))
    # None in underflow
    assert np.all(idx != 0)
    # None in overflow
    assert np.all(idx != len(edges))

    df[f"PNN{mass}Bin"] = idx.astype(np.uint8)

df.to_hdf(args.outfile, "df_hadhad", complevel=9, format="table")
