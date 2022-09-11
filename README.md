
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

This example illustrate how the toa package can be used to identify differentially expressed genes, generate a reference dataset for Transcript Origin Analysis (TOA), and run a (non)bootstrapped TOA.

``` r
library(toa)

#load example data
data("Chang")
data("epith_mesen_ref_raw")

#get differentially expressed genes
DEG_result <- get_DEG(x = Chang$stress, 
                      genes = subset(Chang, select = -stress), 
                      foldThreshDEG = 1.25)
table(DEG_result$DEG)

#get diagnosticity scores
toa_ref_epith_mesen <- get_toa_ref(gene_symbols = epith_mesen_ref_raw[,1],
                                   exp_treatment = epith_mesen_ref_raw[,2:11],
                                   exp_control = epith_mesen_ref_raw[,12:21])

#toa (does not produce bootstrap estimates of mean diagnosticity scores)
toa_result <- toa(x = Chang$stress, 
                  genes = subset(Chang, select = -stress), 
                  toa_ref = toa_ref_epith_mesen, 
                  foldThreshDEG = 1.25)
toa_result$df_results

#toa_boot (produces bootstrap estimates of mean diagnosticity scores)
toa_boot_result <- toa_boot(x = Chang$stress, 
                            genes = subset(Chang, select = -stress), 
                            toa_ref = toa_ref_epith_mesen, 
                            foldThreshDEG = 1.25)
toa_boot_result$df_results
```
