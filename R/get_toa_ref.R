#' @export
get_toa_ref <- function(gene_symbols, exp_treatment, exp_control, logZeroSub = .001, na.rm = TRUE){

  #generate copies of the data
  temp <- cbind(gene_symbols, exp_treatment, exp_control)
  temp2 <- temp

  #log transform all values in one copy of the data (worth noting log values below the logZeroSubstituteExpressionValue are replaced by the logZeroSubstituteExpressionValue)
  temp2[,-1] <- apply(temp2[,-1], c(1,2), function(x) max(log(x), logZeroSub))

  #get column positions for treatment and control data
  treat_cols <- 2:(ncol(exp_treatment) + 1)
  contr_cols <- (max(treat_cols) + 1):(max(treat_cols)  + ncol(exp_control))

  #compute means, sds, and  the diagnosticity scores
  temp$m1 <- apply(temp[,treat_cols], 1, function(x) mean(x, na.rm = na.rm))
  temp$m2 <- apply(temp[,contr_cols], 1, function(x) mean(x, na.rm = na.rm))
  temp$sd1 <- apply(temp[,treat_cols], 1, function(x) sd(x, na.rm = na.rm))
  temp$sd2 <- apply(temp[,contr_cols], 1, function(x) sd(x, na.rm = na.rm))
  temp$ln_m1 <- apply(temp2[,treat_cols], 1, function(x) mean(x, na.rm = na.rm))
  temp$ln_m2 <- apply(temp2[,contr_cols], 1, function(x) mean(x, na.rm = na.rm))
  temp$ln_sd1 <- apply(temp2[,treat_cols], 1, function(x) sd(x, na.rm = na.rm))
  temp$ln_sd2 <- apply(temp2[,contr_cols], 1, function(x) sd(x, na.rm = na.rm))

  temp$DiagnosticityScoresLinear = (temp$m2 - temp$m1) / (temp$sd1^2 + temp$sd2^2)^.5
  temp$DiagnosticityScoresLog = (temp$ln_m2 - temp$ln_m1) / (temp$ln_sd1^2 + temp$ln_sd2^2)^.5

  #aggregate across genes (and fix the name of the gene column after aggregation)
  diagnostic_scores <- aggregate(subset(temp, select = c("DiagnosticityScoresLinear", "DiagnosticityScoresLog")), list(temp[,1]), mean, na.rm = TRUE)
  colnames(diagnostic_scores)[1] <- "gene"

  return(diagnostic_scores)
}
