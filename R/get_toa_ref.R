#' get reference data needed for Transcript Origin Analysis
#'
#' computes diagnosticity scores for a given set of treatment v. control observations, typically replicates.
#'
#' @param gene_symbols a character vector containing gene symbols.
#' @param exp_treatment a numeric matrix of gene expression values, where columns are observations from a treatment/exposure group.
#' @param exp_control a numeric matrix of gene expression values, where columns are observations from a control/reference group.
#' @param logZeroSub the numeric value used to replace expression values that cannot be log2 transformed; .001 by default.
#' @param toupper bool indicating whether gene symbols should be set to uppercase.
#' @return a data frame containing unique gene symbols and average diagnosticity scores.
#' @examples
#' \dontrun{
#' #load example data
#' data("epith_mesen_ref_raw")
#' #get diagnosticity scores
#' toa_ref_epith_mesen <- get_toa_ref(gene_symbols = epith_mesen_ref_raw[,1],
#'                                    exp_treatment = epith_mesen_ref_raw[,2:11],
#'                                    exp_control = epith_mesen_ref_raw[,12:21])
#' }
#' @export
get_toa_ref <- function(gene_symbols, exp_treatment, exp_control, logZeroSub = .001){

  #generate copies of the data
  temp <- cbind(gene_symbols, exp_treatment, exp_control)
  temp2 <- temp

  #log transform all values in one copy of the data (worth noting invalid log values are replaced by the logZeroSub)
  temp2[,-1] <- apply(temp2[,-1], c(1,2), function(x) max(log(x), logZeroSub))

  #get column positions for treatment and control data
  treat_cols <- 2:(ncol(exp_treatment) + 1)
  contr_cols <- (max(treat_cols) + 1):(max(treat_cols)  + ncol(exp_control))

  #compute means, sds, and  the diagnosticity scores
  temp$m1 <- apply(temp[,treat_cols], 1, function(x) mean(x, na.rm = TRUE))
  temp$m2 <- apply(temp[,contr_cols], 1, function(x) mean(x, na.rm = TRUE))
  temp$sd1 <- apply(temp[,treat_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  temp$sd2 <- apply(temp[,contr_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  temp$ln_m1 <- apply(temp2[,treat_cols], 1, function(x) mean(x, na.rm = TRUE))
  temp$ln_m2 <- apply(temp2[,contr_cols], 1, function(x) mean(x, na.rm = TRUE))
  temp$ln_sd1 <- apply(temp2[,treat_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  temp$ln_sd2 <- apply(temp2[,contr_cols], 1, function(x) stats::sd(x, na.rm = TRUE))

  temp$DiagnosticityScoresLinear = (temp$m2 - temp$m1) / (temp$sd1^2 + temp$sd2^2)^.5
  temp$DiagnosticityScoresLog = (temp$ln_m2 - temp$ln_m1) / (temp$ln_sd1^2 + temp$ln_sd2^2)^.5

  #aggregate across genes (and fix the name of the gene column after aggregation)
  diagnostic_scores <- stats::aggregate(subset(temp, select = c("DiagnosticityScoresLinear", "DiagnosticityScoresLog")), list(temp[,1]), mean, na.rm = TRUE)
  colnames(diagnostic_scores)[1] <- "gene"

  #clean gene symbols
  diagnostic_scores$gene <- set_symbols(diagnostic_scores$gene)

  return(diagnostic_scores)
}
