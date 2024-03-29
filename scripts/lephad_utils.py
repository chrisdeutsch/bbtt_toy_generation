import pickle
import numpy as np
import pandas as pd
import uproot

from utils import masspoints


# Bin edges before rebinning (i.e. original histograms)
edgesLephadPreRebin = np.concatenate([
    np.arange(991, dtype=np.float64) / 1000.,
    0.99 + np.arange(1, 101, dtype=np.float64) / 10000.
])


# Bin edges after rebinning
with open("edges_slt.pkl", "rb") as f:
    edgesSLT = pickle.load(f)

with open("edges_ltt.pkl", "rb") as f:
    edgesLTT = pickle.load(f)


def getLephadDf(filename, treename):
    with uproot.open(filename) as f:
        t = f[treename]
        df = t["MVA"].array(library="pd")

    # Convert the 1-level multiindex to regular index
    df.columns = [col for col, in df.columns]

    # Remove columns that are not needed
    cols_keep = []
    cols_keep += [f"PNN_{mass}" for mass in masspoints]
    cols_keep += ["weight", "BTag", "tau0_truth_match",
                  "b0_truth_match", "b1_truth_match",
                  "mcChannelNumber"]

    cols_drop = [col for col in df.columns if col not in set(cols_keep)]
    df.drop(columns=cols_drop, inplace=True)

    # Harmonize PNN name with hadhad
    df.rename(columns={f"PNN_{mass}": f"PNN{mass}" for mass in masspoints},
              inplace=True)

    # Everything is stored as float, transform where applicable
    df = df.astype({
        "BTag": np.int64,
        "mcChannelNumber": np.int64,
        "b0_truth_match": np.int64,
        "b1_truth_match": np.int64,
        "tau0_truth_match": np.int64,
    })

    # Add sample columns
    df["sample"] = treename

    # Apply selections
    sel = (df["BTag"] == 2)
    df = df.loc[sel]
    df.drop(columns="BTag", inplace=True)

    # Sample for which to keep only true taus
    if treename in {"ttbar", "Wtt", "W"}:
        sel = (np.abs(df["tau0_truth_match"]) == 11) \
            | (np.abs(df["tau0_truth_match"]) == 13) \
            | (np.abs(df["tau0_truth_match"]) == 15)
        df = df.loc[sel]

    return df.copy()
