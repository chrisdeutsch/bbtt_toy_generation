#!/usr/bin/env python
from tqdm import tqdm
import argparse
import h5py
import numpy as np
import pandas as pd

from utils import masspoints


parser = argparse.ArgumentParser()
parser.add_argument("dataframe")
parser.add_argument("-o", "--outfile", required=True)
parser.add_argument("--veto-negative-rates", action="store_true")
args = parser.parse_args()


df = pd.read_hdf(args.dataframe)
df["weight"] = df["weight"].astype(np.float64)
df["weightSquared"] = df["weight"]**2

# Drop data
df = df.loc[df["sample"] != "data"]
if "data" in df["sample"].cat.categories:
    df["sample"] = df["sample"].cat.remove_categories(["data"])

# Apply scale factors
zhf_scale = 1.35
ttbar_scale = 0.97

# Z+HF
mask_zhf = \
    (df["sample"] == "Zttbb") | (df["sample"] == "Zttbc") | (df["sample"] == "Zttcc") \
    | (df["sample"] == "Zbb") | (df["sample"] == "Zbc") | (df["sample"] == "Zcc")
df.loc[mask_zhf, "weight"] *= zhf_scale

# ttbar
mask_ttbar = (df["sample"] == "ttbar") | (df["sample"] == "ttbarFakesMC")
df.loc[mask_ttbar, "weight"] *= ttbar_scale

print("Yield after scale factors:")
print(df.groupby("sample")["weight"].agg(Entries="count", Integral="sum"))


# List of all bins for channel
all_bins = []
for mass in sorted(masspoints):
    unique_bins = sorted(df[f"PNN{mass}Bin"].unique())
    for ibin in unique_bins:
        all_bins.append((mass, ibin))

# Matrices for correlation calculation
nbins = len(all_bins)

l1_mat = np.zeros((nbins, nbins), dtype=np.float64)
l2_mat = np.zeros((nbins, nbins), dtype=np.float64)
l3_mat = np.zeros((nbins, nbins), dtype=np.float64)

l1_sumw2_mat = np.zeros((nbins, nbins), dtype=np.float64)
l2_sumw2_mat = np.zeros((nbins, nbins), dtype=np.float64)
l3_sumw2_mat = np.zeros((nbins, nbins), dtype=np.float64)

for i, (mass_i, ibin_i) in tqdm(enumerate(all_bins), total=nbins):
    in_i = df[f"PNN{mass_i}Bin"] == ibin_i

    for j, (mass_j, ibin_j) in enumerate(all_bins):
        in_j = df[f"PNN{mass_j}Bin"] == ibin_j

        l1_mat[i, j] = df.loc[in_i & ~in_j, "weight"].sum() # In bin i but not in j
        l2_mat[i, j] = df.loc[~in_i & in_j, "weight"].sum() # In bin j but not in i
        l3_mat[i, j] = df.loc[in_i & in_j, "weight"].sum()  # In both bin i and j

        l1_sumw2_mat[i, j] = df.loc[in_i & ~in_j, "weightSquared"].sum() # In bin i but not in j
        l2_sumw2_mat[i, j] = df.loc[~in_i & in_j, "weightSquared"].sum() # In bin j but not in i
        l3_sumw2_mat[i, j] = df.loc[in_i & in_j, "weightSquared"].sum()  # In both bin i and j


cnt_neg1 = np.count_nonzero(l1_mat < 0)
cnt_neg2 = np.count_nonzero(l2_mat < 0)
cnt_neg3 = np.count_nonzero(l3_mat < 0)

if args.veto_negative_rates:
    # Set negative rates to 0
    print(f"Setting {cnt_neg1} negative lambda1's to 0...")
    l1_mat[l1_mat < 0] = 0.0
    print(f"Setting {cnt_neg2} negative lambda2's to 0...")
    l2_mat[l2_mat < 0] = 0.0
    print(f"Setting {cnt_neg3} negative lambda3's to 0...")
    l3_mat[l3_mat < 0] = 0.0
else:
    print("Not vetoing negative rates")
    print(f"Negative lambda 1's: {cnt_neg1}")
    print(f"Negative lambda 2's: {cnt_neg2}")
    print(f"Negative lambda 3's: {cnt_neg3}")


corr_mat = l3_mat / np.sqrt((l1_mat + l3_mat) * (l2_mat + l3_mat))

# Save to HDF5
with h5py.File(args.outfile, "w") as fout:
    fout.create_dataset("bin_labels", data=np.array(all_bins, dtype=int))
    fout.create_dataset("corr", data=corr_mat)
    fout.create_dataset("l1", data=l1_mat)
    fout.create_dataset("l2", data=l2_mat)
    fout.create_dataset("l3", data=l3_mat)
