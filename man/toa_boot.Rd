% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toa_boot.R
\name{toa_boot}
\alias{toa_boot}
\title{bootstrapped Transcript Origin Analysis}
\usage{
toa_boot(toa_result, n_boot = 200, verbose = TRUE)
}
\arguments{
\item{toa_result}{a toa result object produced using \code{toa()}.}

\item{n_boot}{the number of bootstrap sample to run.}

\item{verbose}{a bool indicating whether a progress bar should be shown.}
}
\value{
a list object containing:
\enumerate{
\item a results data frame.
\item a matrix containing bootstrapped mean diagnosticity score estimates \code{boot}.
\item a list of bootstrapped statistics.
\item a list of input arguments.
}
}
\description{
performs transcript origin analysis, including bootstrapped estimation of mean
diagnosticity scores; will use n-1 available CPU cores by default.
}
\examples{
\dontrun{
#load example data
data("TAU_Trials3_Gene_CPM_Log2")
data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
data("HumanM1M2_3Reps_Martinez")

#get DEG
DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
                      regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
                      foldThreshDEG = 2)

#get toa
toa_result <- toa(DEG_result = DEG_result,
                  ref = HumanM1M2_3Reps_Martinez,
                  type_1_cols = 2:4,
                  type_2_cols = 5:7)

#get bootstrapped stats
toa_boot_result <- toa_boot(toa_result)
toa_boot_result$df_results
}
}
