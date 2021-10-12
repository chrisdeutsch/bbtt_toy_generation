import numpy as np
import pandas as pd
import pickle
import uproot

from utils import masspoints


# Bin edges before rebinning (i.e. original histograms)
edgesHadhadPreRebin = np.linspace(0, 1, 1001, dtype=np.float64)


# Bin edges after rebinning
with open("edges_hadhad.pkl", "rb") as f:
    edgesHadhad = pickle.load(f)


def getHadhadDf(filename, treename):
    variables = ["weight", "same_sign", "n_btag", "is_sr"]
    variables += ["b0_truth_match", "b1_truth_match", "mc_channel", "run_number"]
    variables += [f"PNN{mass}" for mass in masspoints]

    with uproot.open(filename) as f:
        t = f[treename]
        df = t.arrays(variables, library="pd")

    # Apply selection
    sel = (~df["same_sign"]) & (df["n_btag"] == 2)
    if treename != "Fake":
        sel = sel & df["is_sr"]

    df = df.loc[sel]
    df.drop(columns=["same_sign", "n_btag", "is_sr"], inplace=True)

    # Sample (re-)naming
    name_mapping = {
        "Zee": "Z",
        "Zmumu": "Z",
        "Ztautau": "Ztt",
    }

    if treename in name_mapping:
        df["sample"] = name_mapping[treename]
    else:
        df["sample"] = treename

    return df.copy()
