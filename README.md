# BBMBC: A Beta-Binomial Framework for Evaluating Rare Memory B-Cell Responses in Early-Phase Vaccine Trials

[![CRAN status](https://www.r-pkg.org/badges/version/BBMBC)](https://CRAN.R-project.org/package=BBMBC)

**BBMBC** is an R package providing statistical tools to model overdispersed frequency data for vaccine-specific memory B cells. It enables researchers to calculate statistical power, estimate overdispersion parameters, and perform likelihood ratio tests tailored to the unique distributions often found in immunology clinical trials.

The methodology implemented in this package is motivated by analytical approaches used for HIV vaccine trial data (e.g., Parks et al., 2025). 

> **Note on Simulations:** All R code used in the simulations for our upcoming paper *(DOI link will be provided soon)* is available in this repository.

## Installation

You can install the released version of **BBMBC** from CRAN:

```R
install.packages("BBMBC")
