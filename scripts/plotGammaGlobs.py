#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

import ROOT as R
R.gROOT.SetBatch(True)


parser = argparse.ArgumentParser()
parser.add_argument("toys")
parser.add_argument("asimov")
parser.add_argument("-c", "--channel", choices=["Hadhad", "SLT", "LTT"], required=True)
parser.add_argument("-m", "--mass", default=None, type=int, required=True)
parser.add_argument("--bin", default=-1, type=int)
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


def getTausWS(args):
    channel_dict = {
        "Hadhad": "hh",
        "SLT": "lh_slt",
        "LTT": "lh_ltt"
    }

    hname = f"tau_{channel_dict[args.channel]}_m{args.mass}"

    fin = R.TFile.Open(args.asimov)
    h_tau = fin.Get(hname)

    tau_ws = np.array([h_tau.GetBinContent(ibin)
                       for ibin in range(1, h_tau.GetNbinsX() + 1)])

    fin.Close()

    return tau_ws


def getGlobsBootstrap(args):
    fin = R.TFile.Open(args.toys)
    tree = fin.Get(f"globs_{args.channel.lower()}") # TODO
    globs = np.array([np.array(event.globs) for event in tree])
    fin.Close()

    return globs


tau_ws = getTausWS(args)
globs_bootstrap = getGlobsBootstrap(args)

with np.printoptions(precision=2):
    print(f"From WS:\n{tau_ws}")
    print(f"Bootstrap:\n{globs_bootstrap.mean(axis=0)}")




bin_idx = np.arange(len(tau_ws))[args.bin]



pois = stats.poisson(mu=tau_ws[bin_idx])
a, b = pois.ppf([0.0005, 0.9995])

x = np.arange(a, b)
y = pois.pmf(x)

bins = x + 0.5


fig, ax = plt.subplots()

ax.hist(globs_bootstrap[:,bin_idx], bins=bins, density=True,
        label=f"Poisson-Bootstrap\n($N = {globs_bootstrap.shape[0]}$)")

ax.plot(x, y, "o",
        label="PMF of $\\mathrm{Pois}(\\tau_{cb})$")

ax.set_xlabel("$m_{cb}$")
ax.set_xlim(bins[0], bins[-1])

ax.set_ylabel("Probability Density / Probability Mass")
ax.set_ylim(0, None)


mean = globs_bootstrap[:, bin_idx].mean()
sem = stats.sem(globs_bootstrap[:, bin_idx])


ax.annotate(f"PNN{args.mass}: Bin {bin_idx + 1}\n"
            f"{args.channel}-channel\n"
            f"$\\tau_{{cb}} = {tau_ws[bin_idx]:.2f}$\n"
            f"$\\left< m_{{cb}} \\right> = {mean:.2f} \\pm {sem:.2f}$",
            xy=(0.03, 0.95),
            xycoords="axes fraction",
            va="top", linespacing=1.8)

ax.legend()

if args.outfile is not None:
    fig.savefig(args.outfile)
