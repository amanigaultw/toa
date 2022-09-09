
<!-- README.md is generated from README.Rmd. Please edit that file -->

# toa

<!-- badges: start -->
<!-- badges: end -->

Offers tools to perform transcriptomics analyses, including Transcript Origin Analysis.

## Installation

You can install the development version of toa from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("amanigaultw/toa")
```

## Example

This is a basic example which shows you how to perform transcript origin analysis using a sample gene expression dataset (Chang) and a sample reference dataset (epith_mesen_ref).

``` r
library(toa)

#load example data
data("Chang")
data("epith_mesen_ref")

#get diagnosticity scores
toa_ref <- get_toa_ref(gene_symbols = epith_mesen_ref[,1], 
                       exp_treatment = epith_mesen_ref[,2:11],
                       exp_control = epith_mesen_ref[,12:21])

#get differentially expressed genes as a function of the stress predictor
DEG_result <- get_DEG(x = Chang$stress, 
                      genes = subset(Chang, select = -stress),
                      foldThreshDEG = 1.25)

#get mean diagnosticity scores for these differentially expressed genes
test1 <- toa(x = Chang$stress, 
             genes = subset(Chang, select = -stress), 
             toa_ref = toa_ref,
             foldThreshDEG = 1.25)

test1$df.means
```
