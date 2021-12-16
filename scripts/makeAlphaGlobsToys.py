#!/usr/bin/env python
import argparse
import numpy as np
from scipy import stats
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("workspaces", nargs="+")
parser.add_argument("-o", "--outfile", default="alphas.root")
args = parser.parse_args()

rng = np.random.default_rng(64782119739)

trunc_norm = stats.truncnorm(-5, 5)


import ROOT as R
R.gROOT.SetBatch(True)

def get_globs(filename):
    f = R.TFile.Open(filename)

    w = f.Get("combined")
    model = w.obj("ModelConfig")

    globs = set()
    for param in model.GetGlobalObservables():
        name = param.GetName()
        if not name.startswith("nom_alpha_"):
            continue

        globs.add(name)

    f.Close()
    return globs


all_globs = set()
for fn in args.workspaces:
    all_globs = all_globs | get_globs(fn)


# Filter global observables
# all_globs_filtered = set()
# for glob in sorted(all_globs):
#     if "ATLAS_LUMI" in glob \
#        or "SysEG" in glob \
#        or "SysEL" in glob \
#        or "SysFT" in glob \
#        or "SysJET" in glob \
#        or "SysMET" in glob \
#        or "SysMUON" in glob \
#        or "SysPRW" in glob \
#        or "SysTAUS" in glob:
#         all_globs_filtered.add(glob)

# all_globs = all_globs_filtered

for glob in sorted(all_globs):
    print(glob)


f = R.TFile.Open(args.outfile, "RECREATE")
tree = R.TTree("globs_alphas", "globs_alphas")
tree.SetDirectory(f)

arrays = {}
branches = {}
for glob in sorted(all_globs):
    arrays[glob] = np.zeros(1, dtype=np.float32)
    branches[glob] = tree.Branch(glob, arrays[glob], f"{glob}/F")


for i in tqdm(range(20000)):
    for glob in sorted(all_globs):
        arrays[glob][0] = trunc_norm.rvs(random_state=rng)

    tree.Fill()

tree.Write()
f.Close()
