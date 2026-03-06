# Area-Level Multi-Type Model

This repository contains R code for an area-level Bayesian multi-type model that jointly analyzes a Gaussian outcome and a Binomial outcome using shared spatial structure.

The code is designed for small area estimation settings where related outcomes can borrow strength through a common latent spatial effect. In this repository, the empirical application focuses on county-level analysis using ACS-based variables such as income and poverty rate.

## Repository Structure

- `packages.r`  
  Installs and loads the required R packages.

- `functions.r`  
  Contains the main utility functions and MCMC samplers, including:
  - shared random effect sampler for the joint model
  - univariate Gaussian spatial model
  - univariate Binomial spatial model
  - simulation helpers
  - evaluation functions such as MSE, coverage, and interval score

- `region.r`  
  Runs the empirical area-level analysis and produces maps and variance comparisons.

- `empirical study.r`  
  Runs the simulation study, fits the joint and univariate models, and computes performance metrics.

- `empirical results.r`  
  Summarizes simulation results and generates comparison plots.

## Model Overview

The repository implements a joint Bayesian framework for two types of area-level responses:

- a Gaussian response
- a Binomial response

The two outcomes are linked through a shared spatial random effect. The Binomial component is handled through Pólya–Gamma data augmentation, which allows efficient Gibbs sampling. The code also includes separate univariate models for comparison.

## Requirements

The code is written in R. Required packages are loaded through `packages.r`. Main package dependencies include:

- `Matrix`
- `MASS`
- `mvtnorm`
- `invgamma`
- `BayesLogit`
- `coda`
- `Metrics`
- `sf`
- `spdep`
- `spatialreg`
- `ggplot2`
- `RColorBrewer`
- `tidyverse`
- `dplyr`
- `tidycensus`

## Data

This repository expects external `.RData` files that are not included here. The scripts reference files such as:

- `westnorthc.Rdata`
- `SD cleaned.RData`
- `empirical results.RData`

You will need to prepare these files before running the analysis. They should contain the spatial objects, direct estimates, sampling variances, covariates, and spatial basis matrices used by the scripts.

## How to Run

A typical workflow is:

1. Load required packages:
   ```r
   source("packages.r")
2. Load functions:

source("functions.r")

3. Run the empirical analysis:

source("region.r")

4. Run the simulation study:

source("empirical study.r")

5. Summarize and plot results:

source("empirical results.r")

## Output

The scripts produce:

posterior mean estimates

posterior variance comparisons

spatial maps

MSE comparisons

coverage summaries

interval score summaries

exported figures in .png format

## Notes

The code currently uses Missouri county-level spatial data in the empirical example.

The simulation study compares:

direct estimates

univariate spatial models

the proposed multi-type joint model

The README can be expanded later to include a formal model description, data dictionary, and reproducible example once the input data files are organized for public release.

## Author

Zewei Kong
