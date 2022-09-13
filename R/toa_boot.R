#' bootstrapped Transcript Origin Analysis
#'
#' performs transcript origin analysis, including bootstrapped estimation of mean diagnosticity scores; will use n-1 available CPU cores by default.
#'
#' @param toa a toa result object produced using \code{toa()}
#' @param n_boot the number of bootstrap sample to run
#' @param progress a bool indicating whether a progress bar should be shown
#' @return a list object containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a matrix containing bootstrapped mean diagnosticity score estimates \code{boot}.
#' }
#' @examples
#' #load example data
#' data("Chang")
#' data("epith_mesen_ref_raw")
#' #get diagnosticity scores
#' toa_ref_epith_mesen <- get_toa_ref(gene_symbols = epith_mesen_ref_raw[,1],
#'                                    exp_treatment = epith_mesen_ref_raw[,2:11],
#'                                    exp_control = epith_mesen_ref_raw[,12:21])
#' #toa
#' toa_result <- toa(x = Chang$stress,
#'                   genes = subset(Chang, select = -stress),
#'                   toa_ref = toa_ref_epith_mesen,
#'                   foldThreshDEG = 1.25)
#' #toa_boot
#' #toa_boot_result <- toa_boot(toa = toa_result,
#' #                            foldThreshDEG = 1.25)
#' @export
toa_boot <- function(toa, n_boot = 200, progress = TRUE){

  #verify inputs
  if(class(toa) != "toa" | is.null(toa$df_results) | is.null(toa$df_DEG)) stop("invalid toa parameters; check whether toa() produced a valid result object")

  #re-use toa inputs
  list2env(toa$inputs)

  #instantiate results list
  results <- list(df_results = NULL,
                  boot = NULL,
                  inputs = as.list(environment()))

  #get non bootstrapped means
  non_boot_means <- toa$df_results

  #force cov to matrix if not null
  cov <- cov_to_matrix(cov, x)

  #reset any previous multithreading settings
  env <- utils::getFromNamespace(".foreachGlobals", "foreach")
  rm(list=ls(name=env), pos=env)

  #setup parallel backend to use multiple processors
  cl = parallel::makeCluster(parallel::detectCores()[1]-1) #use all except 1 core
  doParallel::registerDoParallel(cl)
  '%dopar%' <- foreach::'%dopar%'

  #initiate progress bar
  if(progress == TRUE){
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = n_boot, style = 3)
    fprogress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = fprogress)
  }

  #bootstrap loop (multithreaded)
  Vboot = foreach::foreach(i=1:n_boot, .combine='c', .inorder=FALSE, .export=c('get_DEG'), .options.snow = opts) %dopar% {

    #assign default bootstrap results
    bootstrapped_results <- rep(NA, 6)

    #re-sample rows from inputs
    resampled_rows <- sample(seq_len(length(x)), length(x), replace = TRUE)
    resampled_x <- x[resampled_rows]
    resampled_genes <- genes[resampled_rows, ]
    resampled_cov <- cov[resampled_rows, ]

    #get toa means
    temp <- toa_lite(resampled_x, resampled_genes, toa_ref, resampled_cov, foldThreshDEG)

    #update bootstrap results if toa_lite() was successful
    if(!is.null(temp))  bootstrapped_results <- temp$means

    bootstrapped_results
  }

  #stop multithreading
  parallel::stopCluster(cl)

  #close progress bar
  if(progress == TRUE){
    close(pb)
  }

  #create matrix of bootstrapped results (columns are boostrap resamples; rows correspond to variables)
  boot = matrix(Vboot, nrow = nrow(non_boot_means))
  row.names(boot) = row.names(non_boot_means)

  #remove bootstrapped resamples that did not produce valid values
  boot = boot[ , colSums(is.na(boot)) == 0]

  #
  N_valid_boot_resamples = ncol(boot)

  #compute bootstrapped estimates and stats
  boot_mean <- apply(boot[,], 1, function(x) mean(x, na.rm = TRUE))
  boot_se <- apply(boot[,], 1, function(x) stats::sd(x, na.rm = TRUE))
  boot_z = boot_mean / boot_se
  boot_pValue = 2*stats::pnorm(q = abs(boot_z), lower.tail = FALSE)

  #create results dataframe
  df_results <- data.frame(gene_count_total = non_boot_means$total_genes,
                           gene_count_used_for_TOA = non_boot_means$used_genes,
                           non_boot_mean = non_boot_means$means,
                           boot_mean = boot_mean,
                           boot_se = boot_se,
                           boot_z = boot_z,
                           boot_pValue = boot_pValue)

  #print results
  print(paste0(N_valid_boot_resamples, " out of ",  n_boot, " bootstrap resamples produced valid mean diagnosticity scores"))
  print(df_results)

  #if no failure up to this point, update the results list with computed values
  results$df_results <- df_results
  results$boot <- boot

  return(results)
}

