#!/usr/bin/env python
import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("-o", "--outfile", required=True)
parser.add_argument("-m", "--masses", nargs="+", type=int, default=[300, 325])
parser.add_argument("-s", "--scale", type=float, default=1.0)
args = parser.parse_args()


with h5py.File(args.infile, "r") as fin:
    bin_labels = np.array(fin.get("bin_labels"))
    corr = np.array(fin.get("corr"))


def submatrix(masses):
    """Returns a slice of the large correlation matrix"""
    idx = []
    for i, (mass, ibin) in enumerate(bin_labels):
        if mass in masses:
            idx.append(i)
    idx = np.array(idx)
    return corr[idx[:, np.newaxis], idx[np.newaxis, :]]


mat = submatrix(set(args.masses))
labels = [f"({mass},{ibin})" for mass, ibin in bin_labels if mass in set(args.masses)]
figsize = (args.scale * 6.4, args.scale * 4.8)

fig, ax = plt.subplots(figsize=figsize)
sns.heatmap(100 * mat,
            annot=True, fmt="2.0f", annot_kws={"fontsize": 5},
            vmin=-100, vmax=100, center=0,
            xticklabels=labels,
            yticklabels=labels,
            square=True,
            cbar_kws={"label": r"$\rho$ [%]"},
            ax=ax)
ax.tick_params(axis="both", labelsize=7)
ax.set_xlabel("($m_{X}$ / GeV, bin number)")
ax.set_ylabel("($m_{X}$ / GeV, bin number)")
fig.tight_layout()
fig.savefig(args.outfile)
