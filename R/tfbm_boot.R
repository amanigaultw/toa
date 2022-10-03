#' bootstrapped Transcript Origin Analysis for Transcription Factor Binding Motifs
#'
#' performs a bootstrapped variant of transcript origin analysis aimed at determining whether the
#' frequency of Transcription Factor Binding Motif is significantly increased/decreased
#' among differentially expressed genes. Up to n-1 available CPU cores will be used by default.
#' Bootstrapped estimates of the differentially expressed gene set from \code{toa_boot()} are re-used.
#'
#' @param tfbm_result a tfbm result object produced using \code{tfbm()}.
#' @param verbose a bool indicating whether a progress bar should be shown.
#' @return a list object containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a matrix containing bootstrapped tfbm ratios.
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
#' exp_symbol_col = 1,
#' regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
#' reg_id_col = 1,
#' foldThreshDEG = 2,
#' screenSD = 0,
#' verbose = TRUE)
#'
#' #get toa
#' toa_result <- toa(DEG_result = DEG_result,
#' ref = HumanM1M2_3Reps_Martinez,
#' type_1_cols = 2:4,
#' type_2_cols = 5:7)
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
#'
#' #tfbm boot
#' tfbm_boot_result <- tfbm_boot(tfbm_result)
#' View(tfbm_boot_result$df_results)
#' }
#' @export
tfbm_boot <- function(tfbm_result, verbose = TRUE){

  #instantiate results list
  results <- list(df_results = NULL,
                  boot = NULL,
                  inputs = as.list(environment()))
  class(results) <- "tfbm_boot"

  #get boot mat
  boot_mat <- get_boot_mat(tfbm_result$inputs$toa_boot_result$boot_DEG)

  #get df_ratios
  ratio_dfs <- get_ratio_dfs(boot_mat, tfbm_result$inputs$tfbm_ref, verbose)

  #aggregate across matrices
  arr <- array(unlist(ratio_dfs) , c(nrow(ratio_dfs[[1]]), ncol(ratio_dfs[[1]]), length(ratio_dfs)))
  boot <- apply(arr, 1:2, function(x) mean(x, na.rm = TRUE))
  colnames(boot) <- colnames(ratio_dfs[[1]])

  #compute stats
  N_valid_boot_resamples = apply(boot, 2, function(x) sum(!is.na(x) & !is.nan(x)))
  boot_log2_mean <- tfbm_result$df_results$log2_mean_fold_diff
  boot_log2_se <- apply(boot, 2, function(x) stats::sd(x, na.rm = TRUE))
  boot_log2_z = boot_log2_mean / boot_log2_se
  boot_log2_pValue = 2*stats::pnorm(q = abs(boot_log2_z), lower.tail = FALSE)

  #create results dataframe
  df_results <- data.frame(TFBM = colnames(boot),
                           "Mean fold-difference" = round(2^boot_log2_mean,2),
                           "SE fold-difference" = round(2^boot_log2_se,2),
                           "p SE bootstrap" = round(boot_log2_pValue, 4),
                           "n valid ratios" = tfbm_result$df_results$n_valid_ratios,
                           "Mean log ratio" = boot_log2_mean,
                           "SE_log_ratio" = boot_log2_se)

  #
  results$df_results <- df_results
  results$boot <- boot

  return(results)
}
