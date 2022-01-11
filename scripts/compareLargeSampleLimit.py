#!/usr/bin/env python
from scipy import stats
from tqdm import tqdm
import argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np
import uproot

from utils import masspoints


parser = argparse.ArgumentParser()
parser.add_argument("infile_corr")
parser.add_argument("infile_asimov")
parser.add_argument("-c", "--channel", choices=["Hadhad", "SLT", "LTT"], required=True)
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


rng = np.random.default_rng(96667258605)


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


# Scale up expected yields to be in a ~Gaussian regime
min_mu = mu.min()
sf = np.ceil(20. / min_mu)
print(f"Minimum expected rate: {min_mu:.2f}")
print(f"Apply SF: {sf:.1f}")

mu = sf * mu
print(f"Minimum expected rate after scaling: {mu.min():.2f}")

# Diagonalize correlation matrix
eigval, eigvec = np.linalg.eigh(corr)
with np.printoptions(precision=3, suppress=True):
    print(f"\nEigenvalues:\n{eigval}")

eigval[eigval < 1e-12] = 0.0

# Multivariate normal RVS
gaus = []
pois = []

for batch in tqdm(range(50)):
    rnd = stats.norm.rvs(size=(10000, len(eigval)), random_state=rng)
    rnd *= np.sqrt(eigval)
    # Beware: crazy broadcasting
    rnd = (rnd[:, :, np.newaxis] * eigvec.T[np.newaxis]).sum(axis=1)

    # Use Gaussian approximation
    gaus.append(
        (rnd * np.sqrt(mu) + mu).astype(np.float32)
    )

    # Copula approach
    pois.append(
        stats.poisson.ppf(stats.norm.cdf(rnd), mu=mu).astype(np.float32)
    )

gaus = np.concatenate(gaus)
pois = np.concatenate(pois)


# Summary statistics to compare
gaus_mu = gaus.mean(axis=0)
pois_mu = pois.mean(axis=0)

gaus_dmu = (gaus_mu - mu) / mu
pois_dmu = (pois_mu - mu) / mu

print("Average abs. deviation of mu from target:")
print(f"Normal Approx.: {100 * np.abs(gaus_dmu).mean():.3f} %")
print(f"Copula: {100 * np.abs(pois_dmu).mean():.3f} %")

print("Max abs. deviation of mu from target:")
print(f"Normal Approx.: {100 * np.max(np.abs(gaus_dmu)):.3f} %")
print(f"Copula: {100 * np.max(np.abs(pois_dmu)):.3f} %")


gaus_dvar = ((gaus - mu) / mu).var(axis=0)
pois_dvar = ((pois - mu) / mu).var(axis=0)

#gaus_dvar = (gaus_var - mu) / mu
#pois_dvar = (pois_var - mu) / mu



print("Average abs. deviation of variance from target:")
print(f"Normal Approx.: {100 * np.abs(gaus_dvar).mean():.3f} %")
print(f"Copula: {100 * np.abs(pois_dvar).mean():.3f} %")

print("Max abs. deviation of variance from target:")
print(f"Normal Approx.: {100 * np.max(np.abs(gaus_dvar)):.3f} %")
print(f"Copula: {100 * np.max(np.abs(pois_dvar)):.3f} %")


gaus_corr = np.corrcoef(gaus.T)
pois_corr = np.corrcoef(pois.T)

gaus_dcorr = np.abs(gaus_corr - corr)
pois_dcorr = np.abs(pois_corr - corr)

print("Maximum error (Corr. ME):")
print(f"Normal Approx.: {100 * np.max(gaus_dcorr):.2f} %")
print(f"Copula: {100 * np.max(pois_dcorr):.2f} %")

print("Mean absolute difference (Corr. ME)")
print(f"Normal Approx.: {np.mean(100 * np.abs(gaus_dcorr)):.2f} %")
print(f"Copula: {np.mean(100 * np.abs(pois_dcorr)):.2f} %")


import pdb; pdb.set_trace()

# SF 1
#idx = 143
#bin_range = (-0.5, 14.5)
#bin_num = 15

# SF 4
idx = 143
bin_range = (-0.5, 49.5)
bin_num = 50

# SF 50
#idx = 143
#bin_range = (225.5, 375.5)
#bin_num = 150


fig, ax = plt.subplots()
ax.hist(pois[:, idx], range=bin_range, bins=bin_num, histtype="step", label="Copula approach")
ax.hist(gaus[:, idx], range=bin_range, bins=bin_num, histtype="step", label="Normal approx.")

ax.set_xlabel("Pseudo-Data Yield")
ax.set_ylabel("Toy Experiments")

ax.legend()
fig.savefig("plot.pdf")
