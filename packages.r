## Package Loader
## all required packages listed
req_pkgs <- c(
  # Base math / matrix / sampling
  "Matrix", "MASS", "mvtnorm", "invgamma",
  # Bayesian / Polya-Gamma
  "BayesLogit", "coda", "Metrics",
  # Spatial data handling
  "sf", "spdep", "spatialreg",
  # Visualization
  "ggplot2", "RColorBrewer",
  # Data manipulation
  "tidyverse", "dplyr", "tidycensus", "foreign"
)

# check if all required package installed
for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# load all packages
suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(mvtnorm)
  library(invgamma)
  library(BayesLogit)
  library(coda)
  library(Metrics)
  library(spdep)
  library(spatialreg)
  library(RColorBrewer)
  library(tidyverse)
  library(dplyr)
  library(tidycensus)
  library(foreign)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(rlang)
  library(dplyr)
  library(scales)
})

message("All required packages are loaded and ready.")
