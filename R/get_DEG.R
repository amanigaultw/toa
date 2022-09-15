#' get Differentially Expressed Genes
#'
#' performs regressions to determine which gene show differential expression as a function of x; can adjust for covariates using \code{cov}.
#'
#' @param x a numeric vector of values; the predictor variable of interest.
#' @param genes a numeric matrix of gene expression values, where columns are individual genes.
#' @param cov a numeric vector/matrix of values, where columns are individual covariate measures; if NULL, unadjusted regression coefficients (dif) are produced.
#' @param foldThreshDEG the X-fold expression threshold to be exceeded for a given gene to be considered differentially expressed.
#' @param geneScreenSD minimum gene expression SD; genes that show expression variability below threshold are excluded.
#' @return a data frame containing gene symbols, raw regression coefficients, and a DEG variable; DEG values indicate whether a given gene showed expression values beyond threshold (1 = increased; -1 = decreased; 0 = below threshold).
#' @examples
#' #load example data
#' data("Chang")
#' #get differentially expressed genes
#' DEG_result <- get_DEG(x = Chang$stress,
#'                       genes = subset(Chang, select = -stress),
#'                       foldThreshDEG = 1.25)
#' table(DEG_result$DEG)
#' @export
get_DEG <- function(x, genes, cov = NULL, foldThreshDEG = 1.5, geneScreenSD = 0){

  df_DEG <- NULL

  if(!valid_input_x(x) | !valid_input_genes(genes) | !valid_input_cov(cov)){
    warning("invalid inputs; check warnings")
    return(df_DEG)
  }

  if(is.null(cov)){
    X_vars <- matrix(x)
  }else{
    X_vars <- cbind(x, cov)
  }

  #screen out genes with expression variability below threshold
  keep <- apply(genes, 2, function(x) ifelse(stats::sd(x, na.rm = TRUE) > geneScreenSD, TRUE, FALSE))
  genes <- genes[,keep]

  #run fastLmPure regressions for each gene and save only the regression coefficient of the predictor variable
  dif = numeric(ncol(genes))
  XX <- as.matrix(cbind(intercept = 1, X_vars))
  for (i in seq_len(ncol(genes))){
    dif[i] = RcppArmadillo::fastLmPure(XX, genes[,i])$coefficients[2]
  }

  #generate a dataframe containing gene symbols, regression coefficients and DEG
  df_DEG = data.frame(gene = colnames(genes),
                      dif = dif,
                      DEG = sign(dif)*(abs(dif) > log2(foldThreshDEG)))

  return(df_DEG)
}
