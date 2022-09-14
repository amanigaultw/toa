#' Transcript Origin Analysis for Transcription Factor Binding Motifs
#'
#' performs transcript origin analysis for Transcription Factor Binding Motif without bootstrapped estimates of mean diagnosticity scores.
#'
#' @inheritParams get_DEG
#' @param tfbm_ref a list object containing several tfbm reference data frames; the first column of each data frame must contain gene symbols and be named "gene"
#' @param lite a bool indicating whether a stripped down version of the function should be run; FALSE by default.
#' @return tfbm returns an object of class "tfbm"
#'
#' An object of class "tfbm" is a list containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a raw coefficients data frame.
#'   \item a raw ratios data frame.
#'   \item a list of input arguments.
#' }
#' @examples
#' #load example data
#' data("Chang")
#'
#' #load a tfbm database from another repo (because it is >27MB)
#' #Rfssa::load_github_data("https://github.com/amanigaultw/TELiS/blob/main/HumanTransfacTELiS2019.RData")
#'
#' #tfbm
#' #tfbm_result <- tfbm(x <- Chang[,1],
#' #                    genes <- Chang[,-1],
#' #                    tfbm_ref = HumanTransfacTELiS2019,
#' #                    cov = NULL,
#' #                    foldThreshDEG = 1.25)
#' @export
tfbm <- function(x, genes, tfbm_ref, cov = NULL, foldThreshDEG = 1.5, lite = FALSE){

  #log inputs
  if(lite == FALSE) inputs = as.list(environment())

  #instantiate tfbm result object
  results <- list(df_results = NULL,
                  df_DEG = NULL,
                  df_ratios = NULL,
                  inputs = NULL)
  if(lite == FALSE) results$inputs = inputs
  class(results) <- "tfbm"

  #check inputs
  if(!valid_input_x(x) | !valid_input_genes(genes) | !valid_input_cov(cov) | !valid_input_tfbm_ref(tfbm_ref)){
    warning("invalid inputs; check warnings")
    return(results)
  }

  #get differentially expressed genes
  df_DEG = get_DEG(x, genes, cov, foldThreshDEG)

  if(is.null(df_DEG)){
    warning("get_DEG failed")
    return(results)
  }

  #get tfbm ratios (as up-reg over down-reg)
  df_ratios <- get_ratios(df_DEG, tfbm_ref)

  if(is.null(df_ratios)){
    warning("get_ratios failed")
    return(results)
  }

  nValid <- apply(df_ratios, 2, function(x) sum(!is.na(x)))
  log2mean <- apply(df_ratios, 2, function(x) mean(log2(x), na.rm = TRUE))
  mean <- 2^log2mean

  df_results <- data.frame(TFBM = colnames(df_ratios),
                           n_valid_ratios = nValid,
                           mean_fold_dif = mean,
                           log2_mean_fold_diff = log2mean)

  #if no failure up to this point, update the results list with computed values
  results$df_results <- df_results
  if(lite == FALSE){
    results$df_DEG <- df_DEG
    results$df_ratios <- df_ratios
  }

  return(results)
}
