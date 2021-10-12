# Toy Generation for bbtautau

## Requirements

- ROOT 6.22.06
- Python 3 (see `requirements.txt` for additional packages)


## Step 0: Preparation

1. Create a working directory
2. Get the workspaces and workspace building log files
3. Get ntuples for hadhad, lephad SLT + LTT


Before starting your directory tree should look something like this:
```
TODO
```

## Step 1: Discriminant Binning

Parse the WSMaker workspace building log files to get the bin edges
used for the final discriminant.

```bash
parseBins.py workspaces/logs/build_workspace_2HDM_*.txt
```

Three pickle files should be produced that store the bin edges indices
used for rebinning.


## Step 2: Create Asimov Datasets

Produce histograms of the Asimov dataset providing:

1. The background expectation
2. The effective number of MC events (tau) for the Barlow-Beeston method


```bash
mkdir asimov

for mass in 251 260 280 300 325 350 375 400 450 \
                500 550 600 700 800 900 1000 \
                1100 1200 1400 1600; do
    makeAsimov.py workspaces/${mass}.root -m ${mass} -o asimov/asimov_${mass}.root
done

(
    cd asimov
    [[ -e asimov_merged.root ]] && rm asimov_merged.root
    hadd asimov_merged.root asimov_*.root
)
```

The histograms will store the value of the observables `obs_*`,
effective MC statistics `tau_*` for all channels:
`{hh,lh_slt,lh_ltt,zcr}` and for all discriminants `_m*`.


## Step 3: Make Dataframes

Ntuples are converted to pandas dataframes for ease of use during toy
generation. Additionally, the MVA scores are already discretized
according to the binning used in the final fit.


### Step 3.1: Lephad

```bash
mkdir dataframes

makeLephadDataframes.py \
    ntuples/lephad/SLT_Ntuple_{NoDataFakes,DataFakes}.root \
    -c SLT -o dataframes/dataframe_slt.h5
    
makeLephadDataframes.py \
    ntuples/lephad/LTT_Ntuple_{NoDataFakes,DataFakes}.root \
    -c LTT -o dataframes/dataframe_ltt.h5
```


### Step 3.2: Hadhad

```bash
makeHadhadDataframes.py ntuples/hadhad/mva_ntup.root -o dataframes/dataframe_hadhad.h5
```


## Step 4: Expected Correlation Between Bins in Data

Estimates the expected correlation matrix of data yields per bin using
a bivariate Poisson model.

```bash
mkdir correlation_matrices

makeCorr.py dataframes/dataframe_slt.h5 -o correlation_matrices/corr_slt.h5
makeCorr.py dataframes/dataframe_ltt.h5 -o correlation_matrices/corr_ltt.h5
makeCorr.py dataframes/dataframe_hadhad.h5 -o correlation_matrices/corr_hadhad.h5
```

TODO: make plots


## Step 5: Generate Poisson RVS

This generates correlated Poisson random variables for the observables
in the hadhad / lephad SLT / lephad LTT channels:

```bash
mkdir poisson_rvs

generateFromCorr.py correlation_matrices/corr_slt.h5 asimov/asimov_merged.root \
    -c SLT -o poisson_rvs/rvs_slt.h5

generateFromCorr.py correlation_matrices/corr_ltt.h5 asimov/asimov_merged.root \
    -c LTT -o poisson_rvs/rvs_ltt.h5

generateFromCorr.py correlation_matrices/corr_hadhad.h5 asimov/asimov_merged.root \
    -c Hadhad -o poisson_rvs/rvs_hadhad.h5
```

The pseudo-data for the Z-CR is generated from the workspaces in a
later step (Step 8).


## Step 6: Generate Global Observables (Barlow-Beeston)

This will create multiple root-files containing trees with the values
of the global observables.


```bash
mkdir gamma_globs

makeGammaGlobsToys.py dataframes/dataframe_slt.h5 asimov/asimov_merged.root \
    -o gamma_globs -c SLT

makeGammaGlobsToys.py dataframes/dataframe_ltt.h5 asimov/asimov_merged.root \
    -o gamma_globs -c LTT

makeGammaGlobsToys.py dataframes/dataframe_hadhad.h5 asimov/asimov_merged.root \
    -o gamma_globs -c Hadhad
```


## Step 7: Generate Global Observables (Others)

All other (Gaussian-constrained) global observables are treated as
fully correlated. The script takes *all* workspaces so that it can
figure what the names of all NPs are.

```bash
mkdir other_globs

makeAlphaGlobsToys.py \
 workspaces/{251,260,280,300,325,350,375,400,450,500,550,600,700,800,900,1000,1100,1200,1400,1600}.root \
 -o other_globs/alphas.root
```


## Step 8: Generate Z-CR Toys

The Z-CR toys can be directly generated from the workspaces / Asimov dataset:

```bash
mkdir toys_zcr

makeToysZCR.py asimov/asimov_merged.root -o toys_zcr
```

Will produce the pseudo-dataset to be used in the fit (replacing data)
as well as the global observables related to the gamma NPs.


## Step 9: Build Workspace Inputs

This (mostly technical) step translates the pseudo-data that was
generated in the final fit binning to the initial histogram binning so
that it can be used in the workspace building as if it were real data.

**Pseudo-data:**

```bash
mkdir ws_inputs

cp toys_zcr/pseudodata_ZCR.root ws_inputs/

makePseudoDataHists.py poisson_rvs/rvs_slt.h5 -c SLT -o ws_inputs/
makePseudoDataHists.py poisson_rvs/rvs_ltt.h5 -c LTT -o ws_inputs/
makePseudoDataHists.py poisson_rvs/rvs_hadhad.h5 -c Hadhad -o ws_inputs/
```

**Global observables:**

The global observables are stored as trees where the index of the
entry corresponds to the toy experiment. We can just merge the ROOT
files for all global observables:

```bash
# TODO
for mass in ...; do
    hadd ...
done
```


## Step 10: Fit Toys & Evaluation

This is handled in a separate
[repository](https://gitlab.cern.ch/cdeutsch/bbtt_global_significance/).
