#!/usr/bin/env python3
import argparse
import numpy as np
import pickle
import re


parser = argparse.ArgumentParser()
parser.add_argument("infiles", nargs="+")
parser.add_argument("-o", "--outdir", default="")
args = parser.parse_args()


# Regex patterns
pattern_hh = re.compile("""INFO::Category: In category Region_BMin0_incJet1_distPNN(\\d+)_J2_Y2015_DLLOS_T2_SpcTauHH_L0
.*?
nbin \\d+
""", re.DOTALL | re.MULTILINE)

pattern_slt = re.compile("""INFO::Category: In category Region_BMin0_incJet1_dist(\\d+)_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT0_L1
.*?
nbin \\d+
""", re.DOTALL | re.MULTILINE)

pattern_ltt = re.compile("""INFO::Category: In category Region_BMin0_incJet1_dist(\\d+)_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT1_L1
.*?
nbin \\d+
""", re.DOTALL | re.MULTILINE)

pattern_bins = re.compile(r"Bin (\d+) (\d+)")


# Output dictionaries
edges_hh_dict = {}
edges_slt_dict = {}
edges_ltt_dict = {}

for fn in args.infiles:
    with open(fn, "r") as fin:
        log_content = fin.read()

    mass = None

    # hadhad bins
    m = pattern_hh.search(log_content)
    mass = m.group(1)
    bin_idx_hh = pattern_bins.findall(m.group(0))
    bin_idx_hh = [(int(i), int(j)) for i, j in bin_idx_hh]
    bin_idx_hh = [j for i, j in sorted(bin_idx_hh, key=lambda x: x[0], reverse=True)]

    # SLT bins
    m = pattern_slt.search(log_content)
    assert mass == m.group(1)
    bin_idx_slt = pattern_bins.findall(m.group(0))
    bin_idx_slt = [(int(i), int(j)) for i, j in bin_idx_slt]
    bin_idx_slt = [j for i, j in sorted(bin_idx_slt, key=lambda x: x[0], reverse=True)]

    # LTT bins
    m = pattern_ltt.search(log_content)
    assert mass == m.group(1)
    bin_idx_ltt = pattern_bins.findall(m.group(0))
    bin_idx_ltt = [(int(i), int(j)) for i, j in bin_idx_ltt]
    bin_idx_ltt = [j for i, j in sorted(bin_idx_ltt, key=lambda x: x[0], reverse=True)]

    mass = int(mass)

    # Convert to bin edges (after rebinning)
    edges_hh = (np.array(bin_idx_hh, dtype=np.float64) - 1.) / 1000.

    # lephad bins (before rebinning)
    lh_binning = np.concatenate([
        np.arange(991, dtype=np.float64) / 1000.,
        0.99 + np.arange(1, 101, dtype=np.float64) / 10000.
    ])

    # lephad bins (after rebinning)
    bin_idx_slt = np.array(bin_idx_slt)
    edges_slt = lh_binning[bin_idx_slt - 1]

    bin_idx_ltt = np.array(bin_idx_ltt)
    edges_ltt = lh_binning[bin_idx_ltt - 1]


    edges_hh_dict[mass] = edges_hh
    edges_slt_dict[mass] = edges_slt
    edges_ltt_dict[mass] = edges_ltt


with open("edges_hadhad.pkl", "wb") as fout:
    pickle.dump(edges_hh_dict, fout)

with open("edges_slt.pkl", "wb") as fout:
    pickle.dump(edges_slt_dict, fout)

with open("edges_ltt.pkl", "wb") as fout:
    pickle.dump(edges_ltt_dict, fout)
