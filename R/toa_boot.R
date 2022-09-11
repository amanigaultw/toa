#' @export
toa_boot <- function(x, genes, toa_ref, cov = NULL, foldThreshDEG = 1.5, n_boot = 200, show_progress = TRUE){

  #get non bootstrapped means
  argg <- as.list(environment())
  non_boot_means <- do.call(toa, argg)$df.means

  #reset any previous multithreading settings
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)

  #setup parallel backend to use multiple processors
  cores = parallel::detectCores()
  cl = parallel::makeCluster(cores[1]-1) #use all except 1 core
  doParallel::registerDoParallel(cl)

  #initiate progress bar
  if(show_progress == TRUE){
    #suppressMessages(suppressWarnings(require(doSNOW)))
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = Bootstrap_Samples, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }

  '%dopar%' <- foreach::'%dopar%'

  #bootstrap loop (multithreaded)
  Vboot = foreach::foreach(i=1:Bootstrap_Samples, .combine='c', .inorder=FALSE, .export='gene_TOA_function', .options.snow = opts) %dopar% {

    resampled_rows <- sample(seq_len(nrow(Analysis_Dataframe)), nrow(Analysis_Dataframe), replace = TRUE)
    resampled_Analysis_Dataframe <- Analysis_Dataframe[resampled_rows, ]
    argg$Analysis_Dataframe = resampled_Analysis_Dataframe

    temp <- do.call(gene_TOA_function, argg)

    if("df.means" %in% names(temp)){
      bootstrapped_results <- temp$df.means$means
    }else{
      bootstrapped_results <- rep(NA, 6)
    }
  }

  #stop multithreading
  parallel::stopCluster(cl)

  #close progress bar
  if(show_progress == TRUE){
    close(pb)
  }

  #create matrix of bootstrapped results (columns are boostrap resamples; rows correspond to variables)
  boot = matrix(Vboot, nrow = 6)
  row.names(boot) = row.names(non_bootstrapped_results)

  #remove bootstrapped resamples that did not produce valid values
  boot = boot[ , colSums(is.na(boot)) == 0]

  #
  N_valid_boot_resamples = ncol(boot)

  #compute bootstrapped estimates and stats
  boot_mean <- apply(boot[,], 1, function(x) mean(x, na.rm = TRUE))
  boot_se <- apply(boot[,], 1, function(x) sd(x, na.rm = TRUE))
  boot_z = boot_mean / boot_se
  boot_pValue = 2*pnorm(q= abs(boot_z), lower.tail=FALSE)

  #create results dataframe
  results <- data.frame(
    gene_count_total = non_bootstrapped_results$total_genes,
    #gene_count_used_for_TOA = non_bootstrapped_results$used_genes,
    non_boot_mean = non_bootstrapped_results$means,
    boot_mean = boot_mean,
    boot_se = boot_se,
    boot_z = boot_z,
    boot_pValue = boot_pValue
  )

  #print results
  print(paste0(N_valid_boot_resamples, " out of ",  Bootstrap_Samples, " bootstrap resamples produced valid mean diagnosticity scores"))
  print(results)

  #pass results dataframe to the function
  return(list(boot = boot,
              results = results))
}
