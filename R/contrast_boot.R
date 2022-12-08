#' bootstrap test of gene contrast
#'
#' Computes a mean gene contrast score along with bootstrapped tests statistics. A single \code{toa_boot_result} object or the combination of \code{DEG_result} object and \code{n_boot} are acceptable inputs.
#'
#' @param toa_boot_result a toa_boot result object produced using \code{toa_boot()}.
#' @param DEG_result a DEG result object produced using \code{get_DEG()}.
#' @param n_boot the number of bootstrap sample to run.
#' @param contrast a dataframe containing gene symbols (1st column) and contrast loadings (2nd column).
#' @param verbose a bool indicating whether a progress bar should be shown.
#'
#' @return a list object containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a list of input arguments.
#' }
#' @examples
#' \dontrun{
#' #load example data
#' data("TAU_Trials3_Gene_CPM_Log2")
#' data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
#' data("HumanM1M2_3Reps_Martinez")
#' data("CTRAProinflamAntiviralAntibody")
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
#' #METHOD 1 (leverages existing bootstrap resamples)
#' contrast_boot_result_1 <- contrast_boot(toa_boot_result = toa_boot_result,
#'                                         contrast = CTRAProinflamAntiviralAntibody)
#'
#' #METHOD 2 (generates new bootstrap resamples)
#' contrast_boot_result_2 <- contrast_boot(DEG_result = DEG_result,
#'                                         n_boot = 200,
#'                                         contrast = CTRAProinflamAntiviralAntibody)
#'
#' }
#' @export
contrast_boot <- function(toa_boot_result = NULL,
                          DEG_result = NULL,
                          n_boot = NULL,
                          contrast,
                          verbose = TRUE){

  stopifnot("Invalid Inputs; this function accepts either a single toa_boot_result object as input or both a DEG_result object and an n_boot value." = validInputsContrastBoot(toa_boot_result, DEG_result, n_boot))

  #instantiate boot_contrast result object
  results <- list(df_results = NULL,
                  inputs = as.list(environment()))
  class(results) <- "boot_contrast"

  #source data needed based on inputs
  if(!is.null(toa_boot_result)){
    df_DEG_contrast_only <- merge(toa_boot_result$inputs$toa_result$inputs$DEG_result$df_DEG, contrast, by = 1)
    boot_DEG <- toa_boot_result$boot_DEG
  }else{
    df_DEG_contrast_only <- merge(DEG_result$df_DEG, contrast, by = 1)
    boot_DEG <- get_boot_DEG(DEG_result$analysis_data, DEG_result$inputs$regressor_matrix, DEG_result$inputs$foldThreshDEG, n_boot, verbose)
  }

  #handle the possibility of no match
  if(nrow(df_DEG_contrast_only) == 0){
    print("no gene symbols could be matched between the expression_data and contrast matrix provided.")
    return(results)
  }

  #get mean contrast scores
  loading <- df_DEG_contrast_only$dif * df_DEG_contrast_only[,ncol(df_DEG_contrast_only)]
  mean_contrast <- mean(loading, na.rm = TRUE)

  #get boostrapped mean contrast scores
  boot_mat_dif <- get_boot_mat_dif(boot_DEG)

  #compute other stats
  means <- numeric()
  for(i in 2:ncol(boot_mat_dif)){
    df_DEG <- boot_mat_dif[, c(1, i)]
    df_DEG_contrast_only <- merge(df_DEG, contrast, by = 1)
    loading <- as.numeric(df_DEG_contrast_only[,2]) * as.numeric(df_DEG_contrast_only[,3])
    means[i -1] <- mean(loading, na.rm = TRUE)
  }

  se_contrast <- stats::sd(means, na.rm = TRUE)
  z_contrast <- mean_contrast / se_contrast
  p_value_contrast = 2*stats::pnorm(q = abs(z_contrast), lower.tail = FALSE)

  #generate df_result
  df_results <- data.frame(mean_contrast = mean_contrast,
                           se_contrast = se_contrast,
                           z_contrast = z_contrast,
                           p_value_contrast = p_value_contrast)

  #update the result object
  results$df_results <- df_results

  return(results)
}
