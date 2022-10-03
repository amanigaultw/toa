#' Examine relative frequency of Transcription Factor Binding Motifs
#'
#' performs a variant of transcript origin analysis aimed at determining whether the
#' frequency of Transcription Factor Binding Motif is significantly increased/decreased
#' among differentially expressed genes. Bootstrapped estimates of the differentially expressed
#' gene set from \code{toa_boot()} are re-used.
#'
#' @param toa_boot_result a toa_boot result object produced using \code{toa_boot()}.
#' @param tfbm_ref a list object containing several tfbm reference data frames;
#' the first column of each data frame must contain gene symbols and be named "gene"
#' @return tfbm returns an object of class "tfbm"
#'
#' An object of class "tfbm" is a list containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a raw ratios data frame.
#'   \item a list of input arguments.
#' }
#' @examples
#' \dontrun{
#' #load example data
#' data("TAU_Trials3_Gene_CPM_Log2")
#' data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
#' data("HumanM1M2_3Reps_Martinez")
#'
#' #get DEG
#' DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
#'                       regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
#'                       foldThreshDEG = 2)
#'
#' #get toa
#' toa_result <- toa(DEG_result = DEG_result,
#'                   ref = HumanM1M2_3Reps_Martinez,
#'                   type_1_cols = 2:4,
#'                   type_2_cols = 5:7)
#'
#' #get bootstrapped stats
#' toa_boot_result <- toa_boot(toa_result)
#'
#' #load a tfbm database from another repo (because it is >27MB)
#' library(Rfssa)
#' load_github_data("https://github.com/amanigaultw/TELiS/blob/main/HumanTransfacTELiS2019.RData")
#'
#' #tfbm
#' tfbm_result <- tfbm(toa_boot_result, HumanTransfacTELiS2019)
#' View(tfbm_result$df_results)
#' }
#' @export
tfbm <- function(toa_boot_result, tfbm_ref){

  #instantiate tfbm result object
  results <- list(df_results = NULL,
                  df_ratios = NULL,
                  inputs = as.list(environment()))
  class(results) <- "tfbm"

  #get tfbm ratios (as up-reg over down-reg)
  df_ratios <- get_ratios(toa_boot_result$inputs$toa$inputs$DEG$df_DEG, tfbm_ref)

  nValid <- apply(df_ratios, 2, function(x) sum(!is.na(x)))
  log2mean <- apply(df_ratios, 2, function(x) mean(log2(x), na.rm = TRUE))
  mean <- 2^log2mean

  df_results <- data.frame(TFBM = colnames(df_ratios),
                           n_valid_ratios = nValid,
                           mean_fold_dif = mean,
                           log2_mean_fold_diff = log2mean)

  #if no failure up to this point, update the results list with computed values
  results$df_results <- df_results
  results$df_ratios <- df_ratios

  return(results)
}
