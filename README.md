
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

This is a basic example which shows you how to perform transcript origin analysis using a sample gene expression data set (Chang) and a sample reference data set (epith_mesen_ref).

``` r
library(toa)

#load example data
data("Chang")
data("epith_mesen_ref")

#get diagnosticity scores
scores <- compute_diagnosticity_scores(Dataframe = epith_mesen_ref,
                                       gene_col = 1,
                                       Cell_of_interest_col = 2:11,
                                       Reference_cells_col = 12:21,
                                       logZeroSubstituteExpressionValue = .001,
                                       na.rm = TRUE,
                                       gene_to_upper = TRUE)

#examine DEG as a function of the "stress" predictor.
test1 <- gene_TOA_function(Analysis_Dataframe = Chang,
                           predictor_col = "stress",
                           TOA_Reference_Dataframe = scores,
                           covariate_cols = NULL,
                           gene_cols = NULL,
                           ref_gene_col = "gene",
                           foldThreshDEG = 1.25)
test1$df.means

#get bootstrapped TOA estimates
test2 <- bootstrap_gene_TOA_function(Bootstrap_Samples = 200,
                                     Analysis_Dataframe = Chang,
                                     predictor_col = "stress",
                                     TOA_Reference_Dataframe = scores,
                                     covariate_cols = NULL,
                                     gene_cols = NULL,
                                     ref_gene_col = "gene",
                                     foldThreshDEG = 1.25)
```
