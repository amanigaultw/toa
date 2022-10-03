
<!-- README.md is generated from README.Rmd. Please edit that file -->

# toa

<!-- badges: start -->
<!-- badges: end -->

Offers tools to perform transcriptomic analyses, including Transcript Origin Analysis.

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
data("TAU_Trials3_Gene_CPM_Log2")
data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
data("HumanM1M2_3Reps_Martinez")

#get DEG
DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
                      regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
                      foldThreshDEG = 2)
                      
table(DEG_result$df_DEG$DEG)

#get toa
toa_result <- toa(DEG_result = DEG_result,
                  ref = HumanM1M2_3Reps_Martinez,
                  type_1_cols = 2:4,
                  type_2_cols = 5:7)

#get bootstrapped stats
toa_boot_result <- toa_boot(toa_result)

#format bootstrapped toa results and view
pretty_results <- toa_result_format(toa_boot_result)
View(pretty_results)

```

Next, we follow up the toa with an analysis of transcription factor binding motifs

``` r
#load a tfbm database from another repo (because it is >27MB)
library(Rfssa)
load_github_data("https://github.com/amanigaultw/TELiS/blob/main/HumanTransfacTELiS2019.RData")

#tfbm
tfbm_result <- tfbm(toa_boot_result, HumanTransfacTELiS2019)

#tfbm boot
tfbm_boot_result <- tfbm_boot(tfbm_result)
head(tfbm_boot_result$df_results)
```

