geneScreen <- function(expression_data, exp_symbol_col = 1, screenSD = 0, verbose = FALSE){

  #generate a vector indicating whether row SD is below or equal to threshold
  keep <- apply(expression_data[,-exp_symbol_col], 1, function(x) ifelse(stats::sd(x) <= screenSD, 0, 1))

  if(verbose){
    print(paste0(table(keep)[names(table(keep)) == 1], " genes passed geneScreen()"))
  }

  return(expression_data[keep == 1,])
}


pre_DEG_merge <- function(expression_data_screened, regressor_matrix, exp_symbol_col = 1, reg_id_col = 1){

  #transpose and set gene symbols as column names
  analysis_data <- t(expression_data_screened)
  colnames(analysis_data) <- analysis_data[exp_symbol_col,]

  #ensure that using gene symbols as column names does not lead to some R auto-formating issues.
  if(!all.equal(unname(colnames(analysis_data)), unname(analysis_data[1,]))){
    warning("gene symbols were altered when used as column names by R")
    return(NULL)
  }

  #merge with regressor_matrix
  analysis_data <- merge(regressor_matrix, analysis_data, by.x = reg_id_col, by.y = 0, all.x = TRUE)

  #ensure that genes expression values
  analysis_data[,-c(1:ncol(regressor_matrix))] <- apply(analysis_data[,-c(1:ncol(regressor_matrix))], 2, function(x) as.numeric(x))

  return(analysis_data)
}


get_x_cov <- function(analysis_data, regressor_matrix){

  if(ncol(regressor_matrix) > 2){
    x <- analysis_data[,2]
    cov <- analysis_data[,3:ncol(regressor_matrix)]
  }else{
    x <- analysis_data[,2]
    cov <- NULL
  }

  return(list(x = x,
              cov = cov))
}


get_df_DEG <- function(x, genes, cov = NULL, foldThreshDEG = 1.5){

  df_DEG <- NULL

  if(is.null(cov)){
    X_vars <- matrix(x)
  }else{
    X_vars <- cbind(x, cov)
  }

  #run fastLmPure regressions for each gene and save only the regression coefficient of the predictor variable
  dif = numeric(ncol(genes))
  SE = numeric(ncol(genes))
  XX <- as.matrix(cbind(intercept = 1, X_vars))
  for (i in seq_len(ncol(genes))){
    temp = RcppArmadillo::fastLmPure(XX, genes[,i])
    dif[i] = temp$coefficients[2]
    SE[i] = temp$stderr[2]
  }

  tStatistic = dif / SE
  pValue = 2 * stats::pt(q = abs(tStatistic), df = nrow(XX) - ncol(XX), lower.tail = F)

  #generate a dataframe containing gene symbols, regression coefficients and DEG
  df_DEG = data.frame(gene = colnames(genes),
                      dif = dif,
                      SE = SE,
                      tStatistic = tStatistic,
                      pValue = pValue,
                      DEG = sign(dif)*(abs(dif) > log2(foldThreshDEG)))

  df_DEG[sapply(df_DEG, is.nan)] <- NA

  return(df_DEG)
}


get_raw_diag_scores <- function(ref, ref_symbol_col, type_1_cols, type_2_cols, logZeroSub){

  #
  ref_log <- ref
  ref_log[,-ref_symbol_col] <- apply(ref_log[,-ref_symbol_col], c(1,2), function(x) max(log(x), logZeroSub))

  #
  m1 <- apply(ref[,type_1_cols], 1, function(x) mean(x, na.rm = TRUE))
  m2 <- apply(ref[,type_2_cols], 1, function(x) mean(x, na.rm = TRUE))
  sd1 <- apply(ref[,type_1_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  sd2 <- apply(ref[,type_2_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  ref$DiagnosticityScoresLinear = (m2 - m1) / (sd1^2 + sd2^2)^.5

  #
  m1 <- apply(ref_log[,type_1_cols], 1, function(x) mean(x, na.rm = TRUE))
  m2 <- apply(ref_log[,type_2_cols], 1, function(x) mean(x, na.rm = TRUE))
  sd1 <- apply(ref_log[,type_1_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  sd2 <- apply(ref_log[,type_2_cols], 1, function(x) stats::sd(x, na.rm = TRUE))
  ref$DiagnosticityScoresLog = (m2 - m1) / (sd1^2 + sd2^2)^.5

  #
  raw_diag_scores <- ref[,c(names(ref)[ref_symbol_col], "DiagnosticityScoresLinear", "DiagnosticityScoresLog")]

  return(raw_diag_scores)
}


get_final_diag_scores <- function(raw_diag_scores, pop_gene_symbols){

  raw_diag_scores_agg <- stats::aggregate(raw_diag_scores[,-1], list(raw_diag_scores[,1]), mean)
  raw_diag_scores_agg_matched <- raw_diag_scores_agg[raw_diag_scores_agg[,1] %in% pop_gene_symbols, ]

  #
  pop_mean_diag_score_linear <- mean(raw_diag_scores_agg_matched$DiagnosticityScoresLinear)
  pop_mean_diag_score_log <- mean(raw_diag_scores_agg_matched$DiagnosticityScoresLog)

  #computer adjusted diagnosticity scores (adjusting for population mean expression)
  diagnosticity_scores <- data.frame(gene = raw_diag_scores_agg_matched[,1],
                                     diag_score_linear = raw_diag_scores_agg_matched$DiagnosticityScoresLinear - pop_mean_diag_score_linear,
                                     diag_score_log = raw_diag_scores_agg_matched$DiagnosticityScoresLog - pop_mean_diag_score_log)

  return(diagnosticity_scores)
}


get_toa_means <- function(diagnosticity_scores, df_DEG){

  #
  df_final <- merge(diagnosticity_scores, df_DEG, by = 1, all.x = TRUE)

  #
  means <- list(mean_toa_upreg_linear = mean(df_final[df_final$DEG == 1, ]$diag_score_linear),
                mean_toa_downreg_linear = mean(df_final[df_final$DEG == -1, ]$diag_score_linear),
                mean_toa_reg_linear = mean(df_final[df_final$DEG != 0, ]$diag_score_linear),
                mean_toa_upreg_log = mean(df_final[df_final$DEG == 1, ]$diag_score_log),
                mean_toa_downreg_log = mean(df_final[df_final$DEG == -1, ]$diag_score_log),
                mean_toa_reg_log = mean(df_final[df_final$DEG != 0, ]$diag_score_log))

  return(means)
}


get_boot_DEG <- function(analysis_data, regressor_matrix, foldThreshDEG, n_boot, verbose){

  #bootstrap DEG
  #reset any previous multithreading settings
  env <- utils::getFromNamespace(".foreachGlobals", "foreach")
  rm(list=ls(name=env), pos=env)
  cl = parallel::makeCluster(parallel::detectCores()[1]-1) #use all except 1 core
  doParallel::registerDoParallel(cl)
  '%dopar%' <- foreach::'%dopar%'

  #initiate progress bar
  if(verbose == TRUE){
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = n_boot, style = 3)
    fprogress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = fprogress)
  }

  boot_DEG = foreach::foreach(i=1:n_boot, .inorder=FALSE, .options.snow = opts, .export = c("get_df_DEG", "get_x_cov")) %dopar% {
    sampled_rows <- sample(1:nrow(analysis_data), nrow(analysis_data), replace = TRUE)
    DEG <- get_df_DEG(x = get_x_cov(analysis_data, regressor_matrix)$x[sampled_rows],
                      cov = as.data.frame(get_x_cov(analysis_data, regressor_matrix)$cov)[sampled_rows,],
                      genes = analysis_data[sampled_rows,-c(1:ncol(regressor_matrix))],
                      foldThreshDEG = foldThreshDEG)
  }

  #stop multithreading
  parallel::stopCluster(cl)
  if(verbose == TRUE){
    close(pb)
  }

  return(boot_DEG)
}


get_boot_means_matrix <- function(diagnosticity_scores, boot_DEG){

  boot_means_matrix <- matrix(nrow = 200, ncol = 6)
  for(i in 1:length(boot_DEG)){

    #
    df_final <- merge(diagnosticity_scores, boot_DEG[[i]], by = 1, all.x = TRUE)

    #
    boot_means_matrix[i, 1] <- mean(df_final[df_final$DEG == 1, ]$diag_score_linear, na.rm = TRUE)
    boot_means_matrix[i, 2] <- mean(df_final[df_final$DEG == -1, ]$diag_score_linear, na.rm = TRUE)
    boot_means_matrix[i, 3] <- mean(df_final[df_final$DEG != 0, ]$diag_score_linear, na.rm = TRUE)
    boot_means_matrix[i, 4] <- mean(df_final[df_final$DEG == 1, ]$diag_score_log, na.rm = TRUE)
    boot_means_matrix[i, 5] <- mean(df_final[df_final$DEG == -1, ]$diag_score_log, na.rm = TRUE)
    boot_means_matrix[i, 6] <- mean(df_final[df_final$DEG != 0, ]$diag_score_log, na.rm = TRUE)
  }

  return(boot_means_matrix)
}


get_boot_stats <- function(boot_means_matrix, toa_means){

  #declare some variables before assigning values from list2env()
  mean_toa_upreg_linear <- mean_toa_downreg_linear <- mean_toa_reg_linear <- NULL
  mean_toa_upreg_log <- mean_toa_downreg_log <- mean_toa_reg_log <- NULL
  boot_se_toa_upreg_linear <- boot_se_toa_downreg_linear <- boot_se_toa_reg_linear <- NULL
  boot_se_toa_upreg_log <- boot_se_toa_downreg_log <- boot_se_toa_reg_log <- NULL
  z_toa_upreg_linear <- z_toa_downreg_linear <- z_toa_reg_linear <- NULL
  z_toa_upreg_log <- z_toa_downreg_log <- z_toa_reg_log <- NULL

  list2env(toa_means, envir=environment())

  boot_mean <- list(boot_mean_toa_upreg_linear = mean(boot_means_matrix[,1], na.rm = TRUE),
                    boot_mean_toa_downreg_linear = mean(boot_means_matrix[,2], na.rm = TRUE),
                    boot_mean_toa_reg_linear = mean(boot_means_matrix[,3], na.rm = TRUE),
                    boot_mean_toa_upreg_log = mean(boot_means_matrix[,4], na.rm = TRUE),
                    boot_mean_toa_downreg_log = mean(boot_means_matrix[,5], na.rm = TRUE),
                    boot_mean_toa_reg_log = mean(boot_means_matrix[,6], na.rm = TRUE))

  boot_se <- list(boot_se_toa_upreg_linear = stats::sd(boot_means_matrix[,1], na.rm = TRUE),
                  boot_se_toa_downreg_linear = stats::sd(boot_means_matrix[,2], na.rm = TRUE),
                  boot_se_toa_reg_linear = stats::sd(boot_means_matrix[,3], na.rm = TRUE),
                  boot_se_toa_upreg_log = stats::sd(boot_means_matrix[,4], na.rm = TRUE),
                  boot_se_toa_downreg_log = stats::sd(boot_means_matrix[,5], na.rm = TRUE),
                  boot_se_toa_reg_log = stats::sd(boot_means_matrix[,6], na.rm = TRUE))

  list2env(boot_se, envir=environment())

  boot_z <- list(z_toa_upreg_linear = mean_toa_upreg_linear / boot_se_toa_upreg_linear,
                 z_toa_downreg_linear = mean_toa_downreg_linear / boot_se_toa_downreg_linear,
                 z_toa_reg_linear = mean_toa_reg_linear / boot_se_toa_reg_linear,
                 z_toa_upreg_log = mean_toa_upreg_log / boot_se_toa_upreg_log,
                 z_toa_downreg_log = mean_toa_downreg_log / boot_se_toa_downreg_log,
                 z_toa_reg_log = mean_toa_reg_log / boot_se_toa_reg_log)

  list2env(boot_z, envir=environment())

  boot_p_value <- list(p_value_upreg_linear = 2*stats::pnorm(q = abs(z_toa_upreg_linear), lower.tail = FALSE),
                       p_value_downreg_linear = 2*stats::pnorm(q = abs(z_toa_downreg_linear), lower.tail = FALSE),
                       p_value_reg_linear = 2*stats::pnorm(q = abs(z_toa_reg_linear), lower.tail = FALSE),
                       p_value_upreg_log = 2*stats::pnorm(q = abs(z_toa_upreg_log), lower.tail = FALSE),
                       p_value_downreg_log = 2*stats::pnorm(q = abs(z_toa_downreg_log), lower.tail = FALSE),
                       p_value_reg_log = 2*stats::pnorm(q = abs(z_toa_reg_log), lower.tail = FALSE))

  return(list(boot_mean = boot_mean,
              boot_se = boot_se,
              boot_z = boot_z,
              boot_p_value = boot_p_value))
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
    temp <- stats::aggregate(temp[,-(1:ncol(df_DEG))], list(temp$DEG), function(x) mean(x, na.rm = TRUE))
    temp[temp[,] == 0] <- NA
    ratios[[i]] <- temp[temp[,1] == 1, -1] / temp[temp[,1] == -1, -1]
    names(ratios)[i] <- names(tfbm_ref)[i]
  }
  df_ratios <- do.call(rbind, ratios)

  return(df_ratios)
}


get_boot_mat <- function(boot_DEG){

  temp <- list(boot_DEG[[1]][,c("gene")])
  for(i in 1:length(boot_DEG)){
    if(stats::sd(boot_DEG[[i]]$dif, na.rm = TRUE) > 0){
      temp[[i + 1]] <- boot_DEG[[i]]$DEG
    }
  }
  boot_mat <- do.call(cbind, temp)
  colnames(boot_mat) <- NULL

  return(boot_mat)
}


get_ratio_dfs <- function(boot_mat, tfbm_ref, verbose){

  n_boot <- ncol(boot_mat) - 1
  i <- NULL

  #reset any previous multithreading settings
  env <- utils::getFromNamespace(".foreachGlobals", "foreach")
  rm(list=ls(name=env), pos=env)
  available_cores <- parallel::detectCores()[1]-1
  if(available_cores > length(tfbm_ref)){
    available_cores <- length(tfbm_ref)
  }
  cl = parallel::makeCluster(available_cores) #use all except 1 core
  doParallel::registerDoParallel(cl)
  '%dopar%' <- foreach::'%dopar%'

  #initiate progress bar
  if(verbose == TRUE){
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = length(tfbm_ref), style = 3)
    fprogress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = fprogress)
  }

  ratio_dfs = foreach::foreach(i=1:length(tfbm_ref), .inorder=FALSE, .options.snow = opts) %dopar% {
    temp <- merge(boot_mat, tfbm_ref[[i]], by = 1)

    ratios <- list()
    for(j in 1:n_boot){
      temp2 <- sapply(split.data.frame(temp[,-(1:ncol(boot_mat))], temp[,j+1]), colMeans)
      temp2[temp2[,] == 0] <- NA
      ratios[[j]] <- temp2[,"1"] / temp2[,"-1"]
    }
    do.call(rbind, ratios)
  }

  #stop multithreading
  parallel::stopCluster(cl)
  if(verbose == TRUE){
    close(pb)
  }

  return(ratio_dfs)
}

validInputsContrastBoot <- function(toa_boot_result, DEG_result, n_boot){

  if(!is.null(toa_boot_result)){
    return(TRUE)
  }

  if(!is.null(DEG_result) & !is.null(n_boot)){
    return(TRUE)
  }

  return(FALSE)
}

get_boot_mat_dif <- function(boot_DEG){

  temp <- list(boot_DEG[[1]][,c("gene")])
  for(i in 1:length(boot_DEG)){
    if(stats::sd(boot_DEG[[i]]$dif, na.rm = TRUE) > 0){
      temp[[i + 1]] <- boot_DEG[[i]]$dif
    }
  }
  boot_mat <- do.call(cbind, temp)
  colnames(boot_mat) <- NULL

  return(boot_mat)
}
