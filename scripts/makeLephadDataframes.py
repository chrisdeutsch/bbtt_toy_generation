#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

from utils import masspoints
from lephad_utils import getLephadDf, addHeavyFlavourSplit, edgesSLT, edgesLTT


parser = argparse.ArgumentParser()
parser.add_argument("ntuple")
parser.add_argument("ntuple_fake")
parser.add_argument("-o", "--outfile", required=True)
parser.add_argument("-c", "--channel", choices=["SLT", "LTT"], required=True)
args = parser.parse_args()


treenames = [
    "data",
    "ttbar",
    "stops", "stopt", "stopWt",
    "W", "Wtt",
    "Z", "Ztt",
    "DY", "DYtt",
    "ZZ", "WZ", "WW",
    "ggZHtautau", "qqZHtautau", "WHtautau", "VBFHtautau", "ggFHtautau",
    "ttH",
    "ggZHbb", "qqZHbb", "WHbb"
]

dfs = []
for tree in treenames:
    dfs.append(getLephadDf(args.ntuple, tree))

dfs.append(getLephadDf(args.ntuple_fake, "Fake"))

df = pd.concat(dfs)
del dfs

# Adding split of Z+jets in bb, bc, cc, bl, cl, l
addHeavyFlavourSplit(df)
df["sample"] = df["sample"].astype("category")

# Yield tables
print("Yield table")
print(
    df.groupby("sample")["weight"].agg(
        Entries="count",
        Integral="sum")
)

# Adding bin indices
for mass in masspoints:
    # Get bin edges for channel
    edges = None
    if args.channel == "SLT":
        edges = edgesSLT[mass]
    elif args.channel == "LTT":
        edges = edgesLTT[mass]
    else:
        raise RuntimeError(f"Unknown channel: {args.channel}")

    # Discretize MVA scores according to binning
    idx = np.digitize(df[f"PNN{mass}"], bins=edges)

    # Fits uint8
    assert np.all((idx >= 0) & (idx <= 255))
    # None in underflow
    assert np.all(idx != 0)
    # None in overflow
    assert np.all(idx != len(edges))

    df[f"PNN{mass}Bin"] = idx.astype(np.uint8)


df_name = None
if args.channel == "SLT":
    df_name = "df_slt"
elif args.channel == "LTT":
    df_name = "df_ltt"
else:
    raise RuntimeError(f"Unknown channel: {args.channel}")

df.to_hdf(args.outfile, df_name, complevel=9, format="table")
