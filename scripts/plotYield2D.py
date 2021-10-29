#!/usr/bin/env python
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("--bin1", nargs=2, type=int, default=[1000, 3])
parser.add_argument("--bin2", nargs=2, type=int, default=[1100, 3])
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


with h5py.File(args.infile) as f:
    bin_labels = np.array(f["bin_labels"])
    (idx1, ), = np.where(np.all(bin_labels == np.array([args.bin1]), axis=1))
    (idx2, ), = np.where(np.all(bin_labels == np.array([args.bin2]), axis=1))

    yield1 = f["poisson_rvs"][:, idx1]
    yield2 = f["poisson_rvs"][:, idx2]


yield1 = yield1.astype(int)
yield2 = yield2.astype(int)


rho = np.corrcoef(yield1, yield2)[0, 1]


max1 = np.percentile(yield1, 99.9, interpolation="lower")
max2 = np.percentile(yield2, 99.9, interpolation="lower")

fig, ax = plt.subplots()
_, _, _, mesh = ax.hist2d(yield1, yield2,
                       bins=[max1, max2],
                       range=[[-0.5, max1 - 0.5], [-0.5, max2 - 0.5]])

ax.set_xlabel(f"Pseudo-Data Yield in Bin {args.bin1[1]} of PNN{args.bin1[0]}")
ax.set_ylabel(f"Pseudo-Data Yield in Bin {args.bin2[1]} of PNN{args.bin2[0]}")
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label("Toy Experiments")

ax.annotate(f"Pearson's $\\rho = {100*rho:.1f} \\%$", xy=(0.52, 0.9), xycoords="figure fraction")

if args.outfile is not None:
    fig.savefig(args.outfile)
