#!/usr/bin/env python
import argparse
import numpy as np
import os
from scipy import stats
from tqdm import tqdm

from utils import masspoints

import ROOT as R
R.gROOT.SetBatch(True)


parser = argparse.ArgumentParser()
parser.add_argument("asimov")
parser.add_argument("-o", "--outdir", default="")
args = parser.parse_args()


rng = np.random.default_rng(45402781074)

fin = R.TFile.Open(args.asimov)

# Z-CR is independent of the signal mass hypothesis so all histograms
# should be identical (check that this is true)
obs_zcr = None
tau_zcr = None

for mass in masspoints:
    # Use the first point as 'default'
    if not obs_zcr:
        obs_zcr = fin.Get(f"obs_zcr_m{mass}").Clone("obs_zcr")
        obs_zcr.SetDirectory(0)
    if not tau_zcr:
        tau_zcr = fin.Get(f"tau_zcr_m{mass}").Clone("tau_zcr")
        tau_zcr.SetDirectory(0)

    # Check that first point agrees with all others
    obs = fin.Get(f"obs_zcr_m{mass}")
    tau = fin.Get(f"tau_zcr_m{mass}")

    nbins = obs.GetNbinsX()
    assert nbins == tau.GetNbinsX()
    assert nbins == obs_zcr.GetNbinsX()
    assert nbins == tau_zcr.GetNbinsX()

    for i in range(nbins + 2):
        assert abs(obs_zcr.GetBinContent(i) - obs.GetBinContent(i)) < 1e-12
        assert abs(tau_zcr.GetBinContent(i) - tau.GetBinContent(i)) < 1e-12

fin.Close()


nbins = obs_zcr.GetNbinsX()
exp = np.array([obs_zcr.GetBinContent(i) for i in range(nbins + 2)])
tau = np.array([tau_zcr.GetBinContent(i) for i in range(nbins + 2)])

pois_exp = stats.poisson(mu=exp)
pois_tau = stats.poisson(mu=tau)


# Pseudo-data (PD)
fn_out = os.path.join(args.outdir, "pseudodata_ZCR.root")
fout = R.TFile.Open(fn_out, "RECREATE")

for i in tqdm(range(20000)):
    pd = pois_exp.rvs(random_state=rng)

    h_pd = obs_zcr.Clone()
    h_pd.Reset()
    h_pd.SetName(f"PseudoData{i}")
    h_pd.SetTitle(f"PseudoData{i}")

    for ibin in range(nbins + 2):
        h_pd.SetBinContent(ibin, pd[ibin])

    h_pd.Write()

fout.Close()


# Global observables (Barlow-Beeston)
fn_out = os.path.join(args.outdir, "toy_globs_ZCR.root")

f = R.TFile.Open(fn_out, "RECREATE")
tree = R.TTree("globs_ZCR", "globs_ZCR")
tree.SetDirectory(f)

# Care: we don't store under-/ overflow bins (hence: len(tau) - 2)
tree_index = np.zeros(1, dtype=np.int32)
tree_globs = np.zeros(len(tau) - 2, dtype=np.float32)

tree.Branch("index", tree_index, "index/I")
tree.Branch("globs", tree_globs, f"globs[{len(tau) - 2}]/F")

for i in tqdm(range(20000)):
    globs = pois_tau.rvs(random_state=rng)
    tree_index[0] = i
    np.copyto(tree_globs, globs[1:-1])
    tree.Fill()

tree.Write()
fout.Close()
