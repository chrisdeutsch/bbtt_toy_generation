#!/usr/bin/env python
from scipy import stats
from tqdm import tqdm
import argparse
import h5py
import numpy as np
import uproot

from utils import masspoints


parser = argparse.ArgumentParser()
parser.add_argument("infile_corr")
parser.add_argument("infile_asimov")
parser.add_argument("-c", "--channel", choices=["Hadhad", "SLT", "LTT"], required=True)
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


# Use different seed for different channels for reproducibility and
# independent RVS
if args.channel == "Hadhad":
    seed = 1234679123
elif args.channel == "SLT":
    seed = 2384761236
elif args.channel == "LTT":
    seed = 9523609123

rng = np.random.default_rng(seed)


# Read correlation matrix
with h5py.File(args.infile_corr, "r") as fin:
    bin_labels = np.array(fin.get("bin_labels"))
    corr = np.array(fin.get("corr"))


# Read Asimov (for expected background)
asimov_hists = {}
with uproot.open(args.infile_asimov) as fin:
    for mass in masspoints:
        if args.channel == "Hadhad":
            hname = f"obs_hh_m{mass}"
        elif args.channel == "SLT":
            hname = f"obs_lh_slt_m{mass}"
        elif args.channel == "LTT":
            hname = f"obs_lh_ltt_m{mass}"
        else:
            raise RuntimeError("Unknown channel")

        asimov_hists[mass] = fin[hname].to_numpy()


# Build vector of expected rates
mu = []
for mass, ibin in bin_labels:
    hist, edges = asimov_hists[mass]
    # Histograms loaded with uproot do not contain over / underflows
    # but these are counted in bin_labels
    mu.append(hist[ibin - 1])

mu = np.array(mu)


# Diagonalize correlation matrix
eigval, eigvec = np.linalg.eigh(corr)
with np.printoptions(precision=3, suppress=True):
    print(f"Eigenvalues: {eigval}")

eigval[eigval < 1e-12] = 0.0


# Multivariate normal RVS
rvs = []
for batch in tqdm(range(50)):
    rnd = stats.norm.rvs(size=(10000, len(eigval)), random_state=rng)
    rnd *= np.sqrt(eigval)
    # Beware: crazy broadcasting
    rnd = (rnd[:, :, np.newaxis] * eigvec.T[np.newaxis]).sum(axis=1)
    rvs.append(rnd)

rvs = np.concatenate(rvs)
print(rvs.shape)

sample_corr = np.corrcoef(rnd.T)
with np.printoptions(precision=3, suppress=True):
    print("Correlation matrix of RVS:")
    print(sample_corr)

    print("\nMean of RVS:")
    print(rvs.mean(axis=0))

    print("\nStd of RVS:")
    print(rvs.std(axis=0))

# Error
dcorr = np.abs(sample_corr - corr)
print(f"Maximum error: {100 * np.max(dcorr):.2f} %")
print(f"Mean absolute error: {np.mean(100 * np.abs(dcorr)):.2f} %")


# Transform multivariate normal to Poisson
rvs = stats.norm.cdf(rvs)
rvs = stats.poisson.ppf(rvs, mu=mu)


# Check summary statistics to ensure that things worked alright
sample_mu = rvs.mean(axis=0)
sample_var = rvs.var(axis=0)

drel_mu = (sample_mu - mu) / mu
drel_var = (sample_var - mu) / mu

print("\nRelative error on mu:")
print(drel_mu)

print("\nRelative error on variance:")
print(drel_var)

sample_corr = np.corrcoef(rvs.T)
dcorr = np.abs(sample_corr - corr)
print(f"Maximum error: {100 * np.max(dcorr):.2f} %")
print(f"Mean absolute difference: {np.mean(100 * np.abs(dcorr)):.2f} %")

if args.outfile:
    with h5py.File(args.outfile, "w") as fout:
        fout.create_dataset("bin_labels", data=bin_labels)
        fout.create_dataset("poisson_rvs", data=rvs, compression="gzip", compression_opts=9)
