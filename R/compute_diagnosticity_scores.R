# ADD DESCRIPTION
compute_diagnosticity_scores <- function(Dataframe, gene_col, Cell_of_interest_col, Reference_cells_col, logZeroSubstituteExpressionValue = .001, na.rm = TRUE, gene_to_upper = TRUE){

  #make sure that gene symbols are all in upper case characters
  if(gene_to_upper) Dataframe[,gene_col] <- toupper(Dataframe[,gene_col])

  #generate copies of the data
  temp <- Dataframe
  temp2 <- Dataframe

  #log transform all values in one copy of the data (worth noting log values below the logZeroSubstituteExpressionValue are replaced by the logZeroSubstituteExpressionValue)
  temp2[,c(Cell_of_interest_col, Reference_cells_col)] <- apply(temp2[,c(Cell_of_interest_col, Reference_cells_col)], c(1,2), function(x) max(log(x), logZeroSubstituteExpressionValue))

  #compute means, sds, and  the diagnosticity scores
  temp$m1 <- apply(temp[,Cell_of_interest_col], 1, function(x) mean(x, na.rm = na.rm))
  temp$m2 <- apply(temp[,Reference_cells_col], 1, function(x) mean(x, na.rm = na.rm))
  temp$sd1 <- apply(temp[,Cell_of_interest_col], 1, function(x) sd(x, na.rm = na.rm))
  temp$sd2 <- apply(temp[,Reference_cells_col], 1, function(x) sd(x, na.rm = na.rm))
  temp$ln_m1 <- apply(temp2[,Cell_of_interest_col], 1, function(x) mean(x, na.rm = na.rm))
  temp$ln_m2 <- apply(temp2[,Reference_cells_col], 1, function(x) mean(x, na.rm = na.rm))
  temp$ln_sd1 <- apply(temp2[,Cell_of_interest_col], 1, function(x) sd(x, na.rm = na.rm))
  temp$ln_sd2 <- apply(temp2[,Reference_cells_col], 1, function(x) sd(x, na.rm = na.rm))

  temp$DiagnosticityScoresLinear = (temp$m2 - temp$m1) / (temp$sd1^2 + temp$sd2^2)^.5
  temp$DiagnosticityScoresLog = (temp$ln_m2 - temp$ln_m1) / (temp$ln_sd1^2 + temp$ln_sd2^2)^.5

  #aggregate across genes (and fix the name of the gene column after aggregation)
  TOA_diagnostic_scores <- aggregate(subset(temp, select = c("DiagnosticityScoresLinear", "DiagnosticityScoresLog")), list(temp[,gene_col]), mean, na.rm = TRUE)
  colnames(TOA_diagnostic_scores)[1] <- "gene"

  return(TOA_diagnostic_scores)
}
