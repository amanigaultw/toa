toa <- function(x, genes, cov = NULL, toa_ref, foldThreshDEG = 1.5){

  #instantiate results list
  results <- list(df.means = NA,
                  df.DEG = NA,
                  df.DEG.ref.matched = NA,
                  df.DEG.ref.matched.expressed = NA)

  #
  if(!valid_input_x(x) | !valid_input_genes(genes) | !valid_input_cov(cov) | !valid_input_toa_ref(toa_ref)){
    warning("invalid inputs; check warnings")
    return(results)
  }

  #get differentially expressed genes
  df.DEG = get_DEG(x, genes, cov, foldThreshDEG)

  if(is.null(df.DEG)){
    warning("ERROR: get_DEG failed")
    return(results)
  }

  #merge TOA diagnosticity score and genes analyses by using the first column of each dataframe as an ID var
  df.DEG.ref.matched = merge(df.DEG, toa_ref, by = "gene")

  #keep only differentially expressed genes
  df.DEG.ref.matched.expressed <- subset(df.DEG.ref.matched, DEG != 0)

  #compute mean diagnosticity scores (if any up/down regulated genes could be matched up to the reference data)
  if(nrow(df.DEG.ref.matched.expressed) > 0){

    #get DEG counts and diagnosticity score means
    counts <- get_gene_counts(df.DEG, df.DEG.ref.matched)
    means <- get_means(df.DEG.ref.matched, df.DEG.ref.matched.expressed)

    #organize results as a dataframe
    df.means <- data.frame(total_genes = counts$gene_count,
                           used_genes = counts$used_gene_count,
                           means = c(up_reg_linear = means$up_reg_linear,
                                     down_reg_linear = means$down_reg_linear,
                                     reg_linear = means$reg_linear,
                                     up_reg_log = means$up_reg_log,
                                     down_reg_log = means$down_reg_log,
                                     reg_log = means$reg_log))

  }else{
    print("ERROR: could not aggregate diagnosticity score across up and down regulated genes; check how many genes are up/down regulated")
    return(results)
  }

  #if no failure up to this point, update the results list with computed values
  results <- list(df.means = df.means,
                  df.DEG = df.DEG,
                  df.DEG.ref.matched = df.DEG.ref.matched,
                  df.DEG.ref.matched.expressed = df.DEG.ref.matched.expressed)

  return(results)
}
