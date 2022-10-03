#' bootstrapped Transcript Origin Analysis
#'
#' performs transcript origin analysis, including bootstrapped estimation of mean
#' diagnosticity scores; will use n-1 available CPU cores by default.
#'
#' @param toa_result a toa result object produced using \code{toa()}.
#' @param n_boot the number of bootstrap sample to run.
#' @param verbose a bool indicating whether a progress bar should be shown.
#' @return a list object containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a matrix containing bootstrapped mean diagnosticity score estimates \code{boot}.
#'   \item a list of bootstrapped statistics.
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
#' toa_boot_result$df_results
#' }
#' @export
toa_boot <- function(toa_result, n_boot = 200, verbose = TRUE){

  #instantiate toa result object
  results <- list(df_results = NULL,
                  boot_DEG = NULL,
                  boot_means_matrix = NULL,
                  boot_stats = NULL,
                  inputs = as.list(environment()))
  class(results) <- "toa_boot"

  #
  analysis_data <- toa_result$inputs$DEG$analysis_data
  regressor_matrix <- toa_result$inputs$DEG$inputs$regressor_matrix
  foldThreshDEG <- toa_result$inputs$DEG$inputs$foldThreshDEG
  diagnosticity_scores <- toa_result$diagnosticity_scores
  toa_means <- toa_result$toa_means

  #
  boot_DEG <- get_boot_DEG(analysis_data, regressor_matrix, foldThreshDEG, n_boot, verbose)
  boot_means_matrix <- get_boot_means_matrix(diagnosticity_scores, boot_DEG)
  boot_stats <- get_boot_stats(boot_means_matrix, toa_means)

  #generate df_results
  df_results <- data.frame(means = unlist(toa_means),
                           boot_means = unlist(boot_stats$boot_mean),
                           boot_se = unlist(boot_stats$boot_se),
                           boot_z = unlist(boot_stats$boot_z),
                           boot_p_value = unlist(boot_stats$boot_p_value))

  #
  results$df_results <- df_results
  results$boot_DEG <- boot_DEG
  results$boot_means_matrix <- boot_means_matrix
  results$boot_stats <- boot_stats

  return(results)
}

