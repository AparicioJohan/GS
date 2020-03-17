# Genomic Selection <img src="man/figures/logo.png" width="120px" align="right"/>

The GS is a r-package that has been created to run cross validation using BGLR - sommer - ASReml packages. This allows you run CV by command line or with a user interface in shiny.

## Installation

``` r
devtools::install_github("AparicioJohan/GS")
```

``` r
source("https://install-github.me/AparicioJohan/GS")
```

## Run the Shiny app

There's only one exported function in the package and it runs the Shiny app:

``` r
GS::GS_shiny()
```

## Command line

The code below shows how you can run the cross validation.

``` r
# geno <- "imputed_rrBLUP.in"
# samp <- "imputed_rrBLUP_samples.txt"
# phen <- "Phenotypic_Analysis.csv"

# prior <- c("ASReml", "RKHS", "sommer")

# out_table <- crossGP(geno, samp, phen, prior, niter = 2, testporc = 0.3)
```
