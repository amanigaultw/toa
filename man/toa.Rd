% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toa.R
\name{toa}
\alias{toa}
\title{Transcript Origin Analysis}
\usage{
toa(
  DEG_result,
  ref,
  type_1_cols,
  type_2_cols,
  ref_symbol_col = 1,
  logZeroSub = 0.001,
  verbose = TRUE
)
}
\arguments{
\item{DEG_result}{a DEG result object produced using \code{get_DEG()}.}

\item{ref}{a data frame in long format containing gene expression value of the reference data frame.
See \code{data("HumanM1M2_3Reps_Martinez")} for an example of appropriate formatting.}

\item{type_1_cols}{an integer vector indicating the position of replicates for the first "cell type"/"subject type".}

\item{type_2_cols}{an integer vector indicating the position of replicates for the second "cell type"/"subject type".}

\item{ref_symbol_col}{an integer indicating the position of the gene symbol column in the reference data frame;
the first column is used by default.}

\item{logZeroSub}{a numeric value used to substitute expression values in the reference data frame that cannot be
log transformed; .001 is used by default.}

\item{verbose}{a bool indicating whether to print function outputs/info to the console.}
}
\value{
toa returns an object of class "toa"

An object of class "toa" is a list containing:
\enumerate{
\item a results data frame.
\item a vector containing gene symbols that were used to compute mean population diagnosticity scores.
\item a data frame containing unadjusted diagnosticity scores.
\item a data frame containing diagnosticity scores adjusted for mean population diagnosticity scores.
\item a list of toa means.
\item a list of input arguments.
}
}
\description{
performs transcript origin analysis without bootstrapped estimates of mean
diagnosticity scores.
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

toa_result$df_results
}
}
