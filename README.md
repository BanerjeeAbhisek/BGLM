# BGLM

This repository contains the R and C++ implementation of BVSG, together with two runnable analysis scripts:

- `ExampleCode.R`: a small simulated example demonstrating how to run BVSG.
- `prostate.R`: the prostate cancer gene-expression and KEGG pathway analysis.

## Repository contents

- `BVSG.cpp` — core C++ implementation.
- `BVSG.r` — R wrapper and plotting functions.
- `ExampleCode.R` — simulation example.
- `prostate.R` — prostate cancer analysis.
- `dat.csv` — prostate cancer gene-expression data.
- `pathway.csv` — KEGG pathway identifiers.
- `extdata/` — KEGG pathway XML files.
- `myImagePlot.r` — heatmap plotting function.

## Required R packages

Install the required CRAN packages:

    install.packages(c(
      "Rcpp",
      "RcppArmadillo",
      "coda",
      "lattice",
      "sna",
      "mvtnorm",
      "MASS",
      "lars",
      "monomvn"
    ))

Install the required Bioconductor packages:

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(c("KEGGgraph", "impute"))

## Running the example

Set the working directory to the repository folder and source the script:

    setwd("path/to/BGLM")
    source("ExampleCode.R")

The script simulates a small dataset, runs BVSG with one MCMC chain, and prints the estimated regression coefficients.

### Estimated runtime

`ExampleCode.R` takes approximately **3 seconds** on the author's machine.

## Running the prostate cancer analysis

Set the working directory to the repository folder and run:

    setwd("path/to/BGLM")
    source("prostate.R")

The required files—`BVSG.r`, `BVSG.cpp`, `dat.csv`, `pathway.csv`, `myImagePlot.r`, and the `extdata/` directory—must remain in the repository folder.

The analysis fits the BVSG model to multiple KEGG pathways and produces pathway heatmaps and network plots. The plots are saved in:

    mygraph.pdf



### Estimated runtime

A complete run of `prostate.R` takes approximately **12–18 hours** on the author's machine.

Runtime may vary depending on processor speed, available memory, R package versions, and the MCMC settings.

The current scripts use:

    Nburn  <- 10000
    Niter  <- 20000
    Nchain <- 1
    Nthin  <- 1

