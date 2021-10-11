#!/usr/bin/env python3
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("workspace")
parser.add_argument("-m", "--mass", type=int, required=True)
parser.add_argument("-o", "--outfile", required=True)
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)


def fixDataset(dataset):
    fixed_dataset = R.RooDataSet("ds", "", dataset.get(0), "weightVar")
    for i in range(dataset.numEntries()):
        argset = dataset.get(i)
        fixed_dataset.add(argset, dataset.weight())

    return fixed_dataset


def getHistFromAsimov(dataset, observable, cut):
    dataset_reduced = dataset.reduce(R.RooArgSet(observable), cut)
    return dataset_reduced.createHistogram("h", observable)


f = R.TFile.Open(args.workspace)
w = f.Get("combined")
model = w.obj("ModelConfig")

# Observables
obs_hh = w.obj(
    f"obs_x_Region_BMin0_incJet1_distPNN{args.mass}_J2_Y2015_DLLOS_T2_SpcTauHH_L0")
obs_lh_slt = w.obj(
    f"obs_x_Region_BMin0_incJet1_dist{args.mass}_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT0_L1")
obs_lh_ltt = w.obj(
    f"obs_x_Region_BMin0_incJet1_dist{args.mass}_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT1_L1")
obs_zcr = w.obj(
    "obs_x_Region_BMin0_incJet1_Y2015_DZllbbCR_T2_L2_distmLL_J2")

# POI
mu = model.GetParametersOfInterest().first()
mu.setVal(0.0)
mu.setConstant()

# Norm factors
zhf_nf = None
ttbar_nf = None
for param in model.GetNuisanceParameters():
    name = param.GetName()
    if name == "ATLAS_norm_Zhf":
        zhf_nf = param
    elif name == "ATLAS_norm_ttbar":
        ttbar_nf = param

zhf_nf.setVal(1.35)
zhf_nf.setConstant()

ttbar_nf.setVal(0.97)
ttbar_nf.setConstant()


asimov = R.RooStats.AsymptoticCalculator.MakeAsimovData(
    model,
    R.RooArgSet(mu, zhf_nf, ttbar_nf),
    model.GetGlobalObservables())

# Makes the weighting work
asimov = fixDataset(asimov)


# Histograms of Asimov observables
h_obs_hh = getHistFromAsimov(
    asimov, obs_hh,
    f"channelCat == channelCat::Region_BMin0_incJet1_distPNN{args.mass}_J2_Y2015_DLLOS_T2_SpcTauHH_L0"
)

h_obs_lh_slt = getHistFromAsimov(
    asimov, obs_lh_slt,
    f"channelCat == channelCat::Region_BMin0_incJet1_dist{args.mass}_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT0_L1"
)

h_obs_lh_ltt = getHistFromAsimov(
    asimov, obs_lh_ltt,
    f"channelCat == channelCat::Region_BMin0_incJet1_dist{args.mass}_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT1_L1"
)

h_obs_zcr = getHistFromAsimov(
    asimov, obs_zcr,
    f"channelCat == channelCat::Region_BMin0_incJet1_Y2015_DZllbbCR_T2_L2_distmLL_J2"
)

# Asimov global observables
pattern = re.compile(
    r"^nom_gamma_stat_Region_BMin0_incJet1_.*"
    r"(DZllbbCR|SpcTauHH|SpcTauLH_Y2015_LTT0|SpcTauLH_Y2015_LTT1)"
    r".*_bin_(\d+)$"
)

gamma_globs = {}
for param in model.GetGlobalObservables():
    name = param.GetName()
    m = pattern.match(name)
    if not m:
        continue

    region, ibin = m.groups()
    ibin = int(ibin)

    if region == "SpcTauHH":
        region = "hh"
    elif region == "SpcTauLH_Y2015_LTT0":
        region = "lh_slt"
    elif region == "SpcTauLH_Y2015_LTT1":
        region = "lh_ltt"
    elif region == "DZllbbCR":
        region = "zcr"
    else:
        raise RuntimeError("Unknown value encountered for region")

    gamma_globs.setdefault(region, []).append((ibin, param.getVal()))

# Make histograms from gamma observables
glob_hists = {}
for key in gamma_globs:
    hbins = [glob for ibin, glob in sorted(gamma_globs[key], key=lambda x: x[0])]
    nbins = len(hbins)

    if key == "zcr":
        h = R.TH1F(f"tau_{key}", "tau_{key}", nbins, 75, 110)
    else:
        h = R.TH1F(f"tau_{key}", "tau_{key}", nbins, 0, 1)

    for i, content in enumerate(hbins, start=1):
        h.SetBinContent(i, content)
        h.SetBinError(i, R.TMath.Sqrt(content))

    glob_hists[key] = h


def writeRootObject(obj, name, directory):
    obj.SetName(name)
    obj.SetTitle(name)

    directory.cd()
    obj.Write()


fout = R.TFile.Open(args.outfile, "RECREATE")

# Observables
writeRootObject(h_obs_hh, f"obs_hh_m{args.mass}", fout)
writeRootObject(h_obs_lh_slt, f"obs_lh_slt_m{args.mass}", fout)
writeRootObject(h_obs_lh_ltt, f"obs_lh_ltt_m{args.mass}", fout)
writeRootObject(h_obs_zcr, f"obs_zcr_m{args.mass}", fout)

# Global Observables
writeRootObject(glob_hists["hh"], f"tau_hh_m{args.mass}", fout)
writeRootObject(glob_hists["lh_slt"], f"tau_lh_slt_m{args.mass}", fout)
writeRootObject(glob_hists["lh_ltt"], f"tau_lh_ltt_m{args.mass}", fout)
writeRootObject(glob_hists["zcr"], f"tau_zcr_m{args.mass}", fout)

fout.Close()
