gene_TOA_function <- function(Analysis_Dataframe, predictor_col, TOA_Reference_Dataframe, covariate_cols = NULL, gene_cols = NULL, ref_gene_col = "gene", foldThreshDEG = 1.5, ...){

  #by default will assume that all columns in the dataframe are genes
  if(is.null(gene_cols)){
    gene_col <- names(Analysis_Dataframe)[!names(Analysis_Dataframe) %in% c(predictor_col, covariate_cols)]
  }

  #safeguard against feeding in constants as predictors/covariates
  predictor_is_a_constant <- FALSE
  if(is.null(covariate_cols) & sd(Analysis_Dataframe[,predictor_col]) == 0){
    predictor_is_a_constant <- TRUE
  }else if (!is.null(covariate_cols)){
    sds <- apply(Analysis_Dataframe[,c(predictor_col, covariate_cols)], 2, function(x) sd(x, na.rm = TRUE))
    if(any(sds == 0)){
      predictor_is_a_constant <- TRUE
    }
  }
  if(predictor_is_a_constant == TRUE){
    print("ERROR: at least one predictor/covaraite variable is a constant")

    results <- list(df.means = NULL,
                    df.DEG = NULL,
                    df.DEG.ref.matched = NULL,
                    df.DEG.ref.matched.expressed = NULL,
                    cellDiagnosticityPopulationMeanScoresLinear = NULL,
                    cellDiagnosticityPopulationMeanScoresLog = NULL)

    return(results)
  }

  #get genes only dataframe and predictor dataframe
  genes = Analysis_Dataframe[,gene_col]
  X_vars = Analysis_Dataframe[,c(predictor_col, covariate_cols)]

  #initiate vectors to hold non-bootstrapped results
  dif = numeric(ncol(genes))

  #generate predictor matrix
  XX <- as.matrix(cbind(intercept = 1, X_vars))

  #run fastLmPure regressions for each gene and save only the regression coefficient of the predictor variable
  for (i in seq_len(ncol(genes))){
    dif[i] = RcppArmadillo::fastLmPure(XX, genes[,i])$coefficients[2]
  }

  #generate a dataframe containing gene symbols, reg coefficient and DEG
  df.temp = data.frame(
    gene = colnames(genes),
    dif = round(dif, 12),
    DEG = sign(dif)*(abs(dif) > log2(foldThreshDEG)))

  #save counts of up and down regulated genes before accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  up_reg_count <- table(df.temp$DEG)["1"]
  down_reg_count <- table(df.temp$DEG)["-1"]
  reg_count <- down_reg_count + up_reg_count
  gene_count <- rep(c(up_reg_count, down_reg_count, reg_count),2)
  names(gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  #merge TOA diagnosticity score and genes analyses by using the first column of each dataframe as an ID var
  df.temp2 = merge(df.temp, TOA_Reference_Dataframe, by.x = "gene", by.y = ref_gene_col)

  #save counts of up and down regulated genes after accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  used_up_reg_count <- table(df.temp2$DEG)["1"]
  used_down_reg_count <- table(df.temp2$DEG)["-1"]
  used_reg_count <- used_down_reg_count + used_up_reg_count
  used_gene_count <- rep(c(used_up_reg_count, used_down_reg_count, used_reg_count),2)
  names(used_gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  #compute population diagnosticity score
  cellDiagnosticityPopulationMeanScoresLinear <- mean(df.temp2$DiagnosticityScoresLinear)
  cellDiagnosticityPopulationMeanScoresLog <- mean(df.temp2$DiagnosticityScoresLog)

  #keep only differentially expressed genes
  df.temp3 <- subset(df.temp2, DEG != 0)

  #compute mean diagnosticity scores (if any up/down regulated genes could be matched up to the reference data)
  if(nrow(df.temp3) > 0){

    # aggregate by DEG (i.e., flagged as up vs. down regulated)
    df.temp.means <- aggregate(df.temp3[, 4:5], list(df.temp3$DEG), mean)

    #compute results
    up_reg_linear = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
    down_reg_linear = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
    reg_linear = mean(df.temp3$DiagnosticityScoresLinear) - cellDiagnosticityPopulationMeanScoresLinear
    up_reg_log = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
    down_reg_log = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
    reg_log = mean(df.temp3$DiagnosticityScoresLog) - cellDiagnosticityPopulationMeanScoresLog

    #organize results as a dataframe
    df.means <- data.frame(total_genes = gene_count,
                           used_genes = used_gene_count,
                           means = c(up_reg_linear = up_reg_linear,
                                     down_reg_linear = down_reg_linear,
                                     reg_linear = reg_linear,
                                     up_reg_log = up_reg_log,
                                     down_reg_log = down_reg_log,
                                     reg_log = reg_log))

  }else{
    df.means <- NULL
    print("ERROR: could not aggregate diagnosticity score across up and down regulated genes; check how many genes are up/down regulated")
  }

  #use a list to return multiple objects
  results <- list(df.means = df.means,
                  df.DEG = df.temp,
                  df.DEG.ref.matched = df.temp2,
                  df.DEG.ref.matched.expressed = df.temp3,
                  cellDiagnosticityPopulationMeanScoresLinear = cellDiagnosticityPopulationMeanScoresLinear,
                  cellDiagnosticityPopulationMeanScoresLog = cellDiagnosticityPopulationMeanScoresLog)

  return(results)
}
