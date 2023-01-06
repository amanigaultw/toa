#' Analyze
#'
#' performs primary analyses from the toa package in succession so as to identify differentially expression genes, then perform
#' transcript origin analysis (or toa), then perform the transcription factor binding motif variant of toa.
#'
#' @inheritParams get_DEG
#' @inheritParams toa
#' @inheritParams toa_boot
#' @inheritParams tfbm
#' @inheritParams tfbm_boot
#' @return toa returns a list object containing:
#' \enumerate{
#'   \item a DEG result object.
#'   \item a toa result object.
#'   \item a toa_boot result object.
#'   \item a tfbm result object.
#'   \item a tfbm_boot result object.
#' }
#' @examples
#' \dontrun{
#' #
#' library(Rfssa)
#'
#' #load example data
#' data("TAU_Trials3_Gene_CPM_Log2")
#' data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
#' data("HumanM1M2_3Reps_Martinez")
#' load_github_data("https://github.com/amanigaultw/TELiS/blob/main/HumanTransfacTELiS2019.RData")
#'
#' #analyze
#' analyze_result <- analyze(expression_data = TAU_Trials3_Gene_CPM_Log2
#'                           regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
#'                           ref = HumanM1M2_3Reps_Martinez,
#'                           type_1_cols = 2:4,
#'                           type_2_cols = 5:7,
#'                           tfbm_ref = HumanTransfacTELiS2019,
#'                           foldThreshDEG = 2)
#' }
#' @export
analyze <- function(expression_data,
                    regressor_matrix,
                    ref,
                    type_1_cols,
                    type_2_cols,
                    tfbm_ref,
                    exp_symbol_col = 1,
                    reg_id_col = 1,
                    ref_symbol_col = 1,
                    foldThreshDEG = 1.5,
                    DEGfun = NULL,
                    screenSD = 0,
                    logZeroSub = .001,
                    n_boot = 200,
                    verbose = TRUE){

  #
  results <- list(DEG_result = NULL,
                  toa_result = NULL,
                  toa_boot_result = NULL,
                  tfbm_result = NULL,
                  tfbm_boot_result = NULL,
                  inputs = as.list(environment()))

  if(verbose){
    cat(paste("Estimating DEG with fold threshold = ", foldThreshDEG, "\n\n\n"))
  }

  #get DEG
  DEG_result <- get_DEG(expression_data = expression_data,
                        exp_symbol_col = exp_symbol_col,
                        regressor_matrix = regressor_matrix,
                        reg_id_col = reg_id_col,
                        foldThreshDEG = foldThreshDEG,
                        DEGfun = DEGfun,
                        screenSD = screenSD,
                        verbose = verbose)

  if(verbose){
    cat("DEG estimation complete \n\n\n")
    cat("Running toa ...")
  }

  #get toa
  toa_result <- toa(DEG_result = DEG_result,
                    ref = ref,
                    type_1_cols = type_1_cols,
                    type_2_cols = type_2_cols,
                    ref_symbol_col = ref_symbol_col,
                    logZeroSub = logZeroSub,
                    verbose = verbose)

  #get bootstrapped stats
  toa_boot_result <- toa_boot(toa_result = toa_result,
                              n_boot = n_boot,
                              verbose = verbose)

  if(verbose){
    cat("toa complete \n\n\n")
    cat("Running tfbm")
  }

  #tfbm
  tfbm_result <- tfbm(toa_boot_result = toa_boot_result,
                      tfbm_ref = tfbm_ref)

  #tfbm boot
  tfbm_boot_result <- tfbm_boot(tfbm_result = tfbm_result,
                                verbose = verbose)

  if(verbose){
    cat("tfbm complete \n\n\n")
  }

  #
  results$DEG_result <- DEG_result
  results$toa_result <- toa_result
  results$toa_boot_result <- toa_boot_result
  results$tfbm_result <- tfbm_result
  results$tfbm_boot_result <- tfbm_boot_result

  return(results)
}
