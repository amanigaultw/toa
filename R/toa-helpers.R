valid_input_x <- function(x){

  if(sum(is.na(x)) > 0 | stats::sd(x) == 0) {
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

  return(TRUE)
}

valid_input_cov <- function(cov){

  if(!is.null(cov)){
    sds <- apply(as.matrix(cov), 2, function(x) stats::sd(x, na.rm = TRUE))

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

valid_input_tfbm_ref <- function(tfbm_ref){

  for(i in 1:length(tfbm_ref)){
    if(colnames(tfbm_ref[[i]])[1] != "gene" | ncol(tfbm_ref[[i]]) < 1){
      warning("invalid tfbm_ref input")
      return(FALSE)
    }
  }

  return(TRUE)
}

get_gene_counts <- function(df_DEG){

  #save counts of up and down regulated genes before accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  up_reg_count <- table(df_DEG$DEG)["1"]
  down_reg_count <- table(df_DEG$DEG)["-1"]
  reg_count <- down_reg_count + up_reg_count
  gene_count <- rep(c(up_reg_count, down_reg_count, reg_count),2)
  names(gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  #save counts of up and down regulated genes after accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  df_DEG.ref.matched <- subset(df_DEG, df_DEG$ref_matched == 1)
  used_up_reg_count <- table(df_DEG.ref.matched$DEG)["1"]
  used_down_reg_count <- table(df_DEG.ref.matched$DEG)["-1"]
  used_reg_count <- used_down_reg_count + used_up_reg_count
  used_gene_count <- rep(c(used_up_reg_count, used_down_reg_count, used_reg_count),2)
  names(used_gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  return(list(gene_count = gene_count,
              used_gene_count = used_gene_count))
}

get_means <- function(df_DEG){

  #subset df_DEG
  df_DEG.ref.matched <- subset(df_DEG, df_DEG$ref_matched == 1)
  df_DEG.ref.matched.expressed <- subset(df_DEG, df_DEG$ref_matched_expressed == 1)

  #compute population diagnosticity score
  cellDiagnosticityPopulationMeanScoresLinear <- mean(df_DEG.ref.matched$DiagnosticityScoresLinear)
  cellDiagnosticityPopulationMeanScoresLog <- mean(df_DEG.ref.matched$DiagnosticityScoresLog)

  # aggregate by DEG (i.e., flagged as up vs. down regulated)
  df.temp.means <- stats::aggregate(df_DEG.ref.matched.expressed[, c("DiagnosticityScoresLinear", "DiagnosticityScoresLog")], list(df_DEG.ref.matched.expressed$DEG), mean)

  #compute results
  up_reg_linear = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
  down_reg_linear = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
  reg_linear = mean(df_DEG.ref.matched.expressed$DiagnosticityScoresLinear) - cellDiagnosticityPopulationMeanScoresLinear
  up_reg_log = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
  down_reg_log = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
  reg_log = mean(df_DEG.ref.matched.expressed$DiagnosticityScoresLog) - cellDiagnosticityPopulationMeanScoresLog

  return(list(up_reg_linear = up_reg_linear,
              down_reg_linear = down_reg_linear,
              reg_linear = reg_linear,
              up_reg_log = up_reg_log,
              down_reg_log = down_reg_log,
              reg_log = reg_log))
}

get_df_means <- function(df_DEG){

  #get DEG counts and diagnosticity score means
  counts <- get_gene_counts(df_DEG)
  means <- get_means(df_DEG)

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

cov_to_matrix <- function(cov, x){

  if(!is.null(cov)){
    cov <- matrix(cov, nrow = length(x))
  }

  return(cov)
}

get_ratios <- function(df_DEG, tfbm_ref){

  df_ratios <- NULL

  temp <- merge(df_DEG[df_DEG$DEG != 0, ], tfbm_ref[[1]], by = "gene")
  if(nrow(temp) == 0){
    warning("no overlap between genes included in tfbm ref and up/down regulated genes")
    return(df_ratios)
  }

  ratios <- list()
  for(i in 1:length(tfbm_ref)){
    temp <- merge(df_DEG[df_DEG$DEG != 0, ], tfbm_ref[[i]], by = "gene")
    temp <- stats::aggregate(temp[,-(1:3)], list(temp$DEG), function(x) mean(x, na.rm = TRUE))
    temp[temp[,] == 0] <- NA
    ratios[[i]] <- temp[temp[,1] == 1, -1] / temp[temp[,1] == -1, -1]
    names(ratios)[i] <- names(tfbm_ref)[i]
  }
  df_ratios <- do.call(rbind, ratios)

  return(df_ratios)
}

# TELiS_to_big_mem <- function(TELiS){
#
#   bm <- list()
#   for(i in 1:length(TELiS)){
#     temp <- TELiS[[i]]
#     row.names(temp) <- temp[,1]
#     temp <- temp[,-1]
#     bm[[i]] <- bigmemory::as.big.matrix(as.matrix(temp))
#   }
#
#   return(bm)
# }
