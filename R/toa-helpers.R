valid_input_x <- function(x){

  if(sum(is.na(x)) > 0 | sd(x) == 0) {
    warning("invalid x input")
    return(FALSE)
  }

  return(TRUE)
}

valid_input_genes <- function(genes){

  if(sum(is.na(genes)) > 0){
    warning("invalid gene input")
    return(FALSE)
  }

  # sds <- apply(as.matrix(genes), 2, function(x) sd(x, na.rm = TRUE))
  # if(any(sds) == 0){
  #   print("invalid gene input")
  #   return(FALSE)
  # }

  return(TRUE)
}

valid_input_cov <- function(cov){

  if(!is.null(cov)){
    sds <- apply(as.matrix(cov), 2, function(x) sd(x, na.rm = TRUE))

    if(sum(is.na(cov)) > 0 | any(sds == 0)){
      warning("invalid cov input")
      return(FALSE)
    }
  }

  return(TRUE)
}

valid_input_toa_ref <- function(toa_ref){

  if(!all.equal(colnames(toa_ref), c("gene", "DiagnosticityScoresLinear", "DiagnosticityScoresLog")) | !nrow(toa_ref) > 0) {
    warning("invalid toa_ref input")
    return(FALSE)
  }

  return(TRUE)
}

get_gene_counts <- function(df.DEG){

  #save counts of up and down regulated genes before accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  up_reg_count <- table(df.DEG$DEG)["1"]
  down_reg_count <- table(df.DEG$DEG)["-1"]
  reg_count <- down_reg_count + up_reg_count
  gene_count <- rep(c(up_reg_count, down_reg_count, reg_count),2)
  names(gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  #save counts of up and down regulated genes after accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  df.DEG.ref.matched <- subset(df.DEG, df.DEG$ref_matched == 1)
  used_up_reg_count <- table(df.DEG.ref.matched$DEG)["1"]
  used_down_reg_count <- table(df.DEG.ref.matched$DEG)["-1"]
  used_reg_count <- used_down_reg_count + used_up_reg_count
  used_gene_count <- rep(c(used_up_reg_count, used_down_reg_count, used_reg_count),2)
  names(used_gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  return(list(gene_count = gene_count,
              used_gene_count = used_gene_count))
}

get_means <- function(df.DEG){

  #subset df.DEG
  df.DEG.ref.matched <- subset(df.DEG, df.DEG$ref_matched == 1)
  df.DEG.ref.matched.expressed <- subset(df.DEG, df.DEG$ref_matched_expressed == 1)

  #compute population diagnosticity score
  cellDiagnosticityPopulationMeanScoresLinear <- mean(df.DEG.ref.matched$DiagnosticityScoresLinear)
  cellDiagnosticityPopulationMeanScoresLog <- mean(df.DEG.ref.matched$DiagnosticityScoresLog)

  # aggregate by DEG (i.e., flagged as up vs. down regulated)
  df.temp.means <- aggregate(df.DEG.ref.matched.expressed[, 4:5], list(df.DEG.ref.matched.expressed$DEG), mean)

  #compute results
  up_reg_linear = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
  down_reg_linear = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
  reg_linear = mean(df.DEG.ref.matched.expressed$DiagnosticityScoresLinear) - cellDiagnosticityPopulationMeanScoresLinear
  up_reg_log = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
  down_reg_log = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
  reg_log = mean(df.DEG.ref.matched.expressed$DiagnosticityScoresLog) - cellDiagnosticityPopulationMeanScoresLog

  return(list(up_reg_linear = up_reg_linear,
              down_reg_linear = down_reg_linear,
              reg_linear = reg_linear,
              up_reg_log = up_reg_log,
              down_reg_log = down_reg_log,
              reg_log = reg_log))
}

get_df_means <- function(df.DEG){

  #get DEG counts and diagnosticity score means
  counts <- get_gene_counts(df.DEG)
  means <- get_means(df.DEG)

  #organize results as a dataframe
  df.means <- data.frame(total_genes = counts$gene_count,
                         used_genes = counts$used_gene_count,
                         means = c(up_reg_linear = means$up_reg_linear,
                                   down_reg_linear = means$down_reg_linear,
                                   reg_linear = means$reg_linear,
                                   up_reg_log = means$up_reg_log,
                                   down_reg_log = means$down_reg_log,
                                   reg_log = means$reg_log))

  return(df.means)
}
