% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfbm_boot.R
\name{tfbm_boot}
\alias{tfbm_boot}
\title{bootstrapped Transcript Origin Analysis for Transcription Factor Binding Motifs}
\usage{
tfbm_boot(tfbm_result, verbose = TRUE)
}
\arguments{
\item{tfbm_result}{a tfbm result object produced using \code{tfbm()}.}

\item{verbose}{a bool indicating whether a progress bar should be shown.}
}
\value{
a list object containing:
\enumerate{
\item a results data frame.
\item a matrix containing bootstrapped tfbm ratios.
\item a list of input arguments.
}
}
\description{
performs a bootstrapped variant of transcript origin analysis aimed at determining whether the
frequency of Transcription Factor Binding Motif is significantly increased/decreased
among differentially expressed genes. Up to n-1 available CPU cores will be used by default.
Bootstrapped estimates of the differentially expressed gene set from \code{toa_boot()} are re-used.
}
\examples{
\dontrun{
#load example data
data("TAU_Trials3_Gene_CPM_Log2")
data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
data("HumanM1M2_3Reps_Martinez")

#get DEG
DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
exp_symbol_col = 1,
regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
reg_id_col = 1,
foldThreshDEG = 2,
screenSD = 0,
verbose = TRUE)

#get toa
toa_result <- toa(DEG_result = DEG_result,
ref = HumanM1M2_3Reps_Martinez,
type_1_cols = 2:4,
type_2_cols = 5:7)

#get bootstrapped stats
toa_boot_result <- toa_boot(toa_result)

#load a tfbm database from another repo (because it is >27MB)
library(Rfssa)
load_github_data("https://github.com/amanigaultw/TELiS/blob/main/HumanTransfacTELiS2019.RData")

#tfbm
tfbm_result <- tfbm(toa_boot_result, HumanTransfacTELiS2019)

#tfbm boot
tfbm_boot_result <- tfbm_boot(tfbm_result)
View(tfbm_boot_result$df_results)
}
}
