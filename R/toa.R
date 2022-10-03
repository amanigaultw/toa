#' Transcript Origin Analysis
#'
#' performs transcript origin analysis without bootstrapped estimates of mean
#' diagnosticity scores.
#'
#' @param DEG_result
#' @param ref a data frame containing gene symbols and diagnosticity scores.
#' @param type_1_cols
#' @param type_2_cols
#' @param ref_symbol_col
#' @param logZeroSub
#' @param verbose
#' @return toa returns an object of class "toa"
#'
#' An object of class "toa" is a list containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a vector containing gene symbols that were used to compute mean population diagnosticity scores.
#'   \item a data frame containing unadjusted diagnosticity scores.
#'   \item a data frame containing diagnosticity scores adjusted for mean population diagnosticity scores.
#'   \item a list of toa means.
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
#' toa_result$df_results
#' }
#' @export
toa <- function(DEG_result, ref, type_1_cols, type_2_cols, ref_symbol_col = 1, logZeroSub = .001, verbose = TRUE){

  #instantiate toa result object
  results <- list(df_results = NULL,
                  pop_gene_symbols = NULL,
                  raw_diag_scores = NULL,
                  diagnosticity_scores = NULL,
                  toa_means = NULL,
                  inputs = as.list(environment()))
  class(results) <- "toa"

  #get DEG
  df_DEG <- DEG_result$df_DEG
  analysis_data <- DEG_result$analysis_data

  #reformat gene symbols for the reference data
  ref[,ref_symbol_col] <- trimws(toupper(ref[,ref_symbol_col]))

  #get pop_gene_symbols
  df_DEG$matched <- ifelse(df_DEG$gene %in% ref[,ref_symbol_col], 1, 0)
  pop_gene_symbols <- df_DEG$gene[df_DEG$matched == 1]
  if(verbose){
    print(paste0(table(df_DEG$matched)[names(table(df_DEG$matched)) == 1], " genes contributed to population mean diagnosticity score computation"))
  }

  #get raw diagnosticity scores
  raw_diag_scores <- get_raw_diag_scores(ref, ref_symbol_col, type_1_cols, type_2_cols, logZeroSub)

  #get population adjusted diagnosticity scores
  diagnosticity_scores <- get_final_diag_scores(raw_diag_scores, pop_gene_symbols)

  #get toa means
  toa_means <- get_toa_means(diagnosticity_scores, df_DEG)

  #generate df_results
  df_results <- data.frame(means = unlist(toa_means))

  #
  results$df_results <- df_results
  results$pop_gene_symbols <- pop_gene_symbols
  results$raw_diag_scores <- raw_diag_scores
  results$diagnosticity_scores <- diagnosticity_scores
  results$toa_means <- toa_means

  #
  return(results)
}
