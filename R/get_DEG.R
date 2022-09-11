get_DEG <- function(x, genes, cov = NULL, foldThreshDEG = 1.5){

  if(!valid_input_x(x) | !valid_input_genes(genes) | !valid_input_cov(cov)) stop("invalid inputs; check warnings")

  if(is.null(cov)){
    X_vars <- matrix(x)
  }else{
    X_vars <- cbind(x, cov)
  }

  #run fastLmPure regressions for each gene and save only the regression coefficient of the predictor variable
  dif = numeric(ncol(genes))
  XX <- as.matrix(cbind(intercept = 1, X_vars))
  for (i in seq_len(ncol(genes))){
    dif[i] = RcppArmadillo::fastLmPure(XX, genes[,i])$coefficients[2]
  }

  #generate a dataframe containing gene symbols, regression coefficients and DEG
  df.DEG = data.frame(
    gene = colnames(genes),
    dif = round(dif, 12),
    DEG = sign(dif)*(abs(dif) > log2(foldThreshDEG)))

  return(df.DEG)
}
