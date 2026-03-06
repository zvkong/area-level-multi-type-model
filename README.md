# Area-Level Multi-Type Model

This repository contains R code for an area-level Bayesian multi-type model for jointly analyzing a Gaussian response and a Binomial response with shared spatial structure.

The main goal of this project is to improve small area estimation by borrowing strength across response types and across neighboring areas. The code includes the proposed joint model, corresponding univariate comparison models, simulation studies, and an empirical application.

## Repository Structure

### `packages.r`
Installs and loads the R packages required for the project.

### `functions.r`
Defines the main functions used throughout the project, including:

- MCMC samplers for the proposed joint model
- Univariate Gaussian spatial model
- Univariate Binomial spatial model
- Simulation and evaluation utilities
- Performance metrics such as mean squared error, coverage, and interval score

### `region.r`
Runs the empirical area-level analysis and generates maps and variance comparisons.

### `empirical study.r`
Runs the simulation study, fits the competing models, and computes performance measures.

### `empirical results.r`
Summarizes simulation results and produces plots for model comparison.

## Model Description

This repository implements a Bayesian framework for jointly modeling two types of area-level outcomes:

- a Gaussian response
- a Binomial response

The two responses are connected through a shared spatial random effect, which allows the model to borrow information across outcomes. The Binomial component is handled using Pólya–Gamma data augmentation, which leads to efficient Gibbs sampling.

In addition to the proposed joint model, the repository also includes separate univariate models for comparison.

## Requirements

The code is written in R.

Please run:

```r
source("packages.r")
````

to install and load the required packages.

The main dependencies include:

* `Matrix`
* `MASS`
* `mvtnorm`
* `invgamma`
* `BayesLogit`
* `coda`
* `Metrics`
* `sf`
* `spdep`
* `spatialreg`
* `ggplot2`
* `RColorBrewer`
* `tidyverse`
* `dplyr`
* `tidycensus`

## Data

Some scripts in this repository rely on external `.RData` files that are not currently included in the repository. These include files such as:

* `westnorthc.Rdata`
* `SD cleaned.RData`
* `empirical results.RData`

These files should contain the spatial objects, direct estimates, variances, covariates, and other intermediate objects required for the empirical study and result summaries.

Because these data files are not included, the repository may not run fully out of the box without additional data preparation.

## Workflow

A typical workflow is:

### 1. Load required packages

```r
source("packages.r")
```

### 2. Load project functions

```r
source("functions.r")
```

### 3. Run the empirical study

```r
source("region.r")
```

### 4. Run the simulation study

```r
source("empirical study.r")
```

### 5. Summarize and visualize results

```r
source("empirical results.r")
```

## Output

The scripts can produce:

* posterior mean estimates
* posterior variance comparisons
* spatial maps
* model comparison plots
* mean squared error summaries
* coverage summaries
* interval score summaries

Some figures are exported as `.png` files.

## Notes

* The repository focuses on area-level joint modeling for mixed response types.
* The empirical analysis is based on county-level spatial data.
* The simulation study compares the proposed joint model with direct estimators and univariate spatial models.
* The repository is currently organized as a research codebase rather than a fully packaged software library.

## Author

Zewei Kong
