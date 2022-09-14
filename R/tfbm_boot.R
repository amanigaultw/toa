#' bootstrapped Transcript Origin Analysis for Transcription Factor Binding Motifs
#'
#' performs transcript origin analysis for Transcription Factor Binding Motif, including bootstrapped estimation of tfbm ratios; will use n-1 available CPU cores by default.
#'
#' @param tfbm a tfbm result object produced using \code{tfbm()}.
#' @param n_boot the number of bootstrap sample to run.
#' @param progress a bool indicating whether a progress bar should be shown.
#' @return a list object containing:
#' \enumerate{
#'   \item a results data frame.
#'   \item a matrix containing bootstrapped mean diagnosticity score estimates \code{boot}.
#'   \item a list of input arguments.
#' }
#' @examples
#' #load example data
#' data("Chang")
#'
#' #load a tfbm database from another repo (because it is >27MB)
#' #Rfssa::load_github_data("https://github.com/amanigaultw/TELiS/blob/main/HumanTransfacTELiS2019.RData")
#'
#' #tfbm
#' #tfbm_result <- tfbm(x <- Chang[,1],
#' #                    genes <- Chang[,-1],
#' #                    tfbm_ref = HumanTransfacTELiS2019,
#' #                    cov = NULL,
#' #                    foldThreshDEG = 1.25)
#' #tfbm_boot
#' #tfbm_boot_results <- tfbm_boot(tfbm_result)
#' @export
tfbm_boot <- function(tfbm, n_boot = 200, progress = TRUE){

  #verify inputs
  if(!inherits(tfbm, "tfbm") | is.null(tfbm$df_results) | is.null(tfbm$df_DEG)) stop("invalid tfbm input parameters; check whether tfbm() produced a valid result object")

  #instantiate results list
  results <- list(df_results = NULL,
                  boot = NULL,
                  inputs = as.list(environment()))

  #get non bootstrapped means
  non_boot_means <- tfbm$df_results

  #reuse tfbm inputs
  x <- genes <- cov <- tfbm_ref <- foldThreshDEG <- NULL #first bind them locally
  list2env(tfbm$inputs, envir = environment())

  #deal with single cov input
  cov <- cov_to_matrix(tfbm$inputs$cov, tfbm$inputs$x)

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
  Vboot = foreach::foreach(i=1:n_boot, .combine='c', .inorder=FALSE, .options.snow = opts) %dopar% {

    #assign default bootstrap results
    bootstrapped_results <- NULL

    #re-sample rows from inputs
    resampled_rows <- sample(seq_len(length(x)), length(x), replace = TRUE)
    resampled_x <- x[resampled_rows]
    resampled_genes <- genes[resampled_rows, ]
    resampled_cov <- cov[resampled_rows, ]

    #get tfbm means
    temp <- toa::tfbm(resampled_x, resampled_genes, tfbm_ref, resampled_cov, foldThreshDEG, TRUE)

    #update bootstrap results if tfbm() was successful
    if(!is.null(temp$df_results))  bootstrapped_results <- temp$df_results$log2_mean_fold_diff

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

  #
  N_valid_boot_resamples = apply(boot, 1, function(x) sum(!is.na(x) & !is.nan(x)))

  #compute bootstrapped estimates and stats
  boot_log2_mean <- apply(boot, 1, function(x) mean(x, na.rm = TRUE))
  boot_log2_se <- apply(boot[,], 1, function(x) stats::sd(x, na.rm = TRUE))
  boot_log2_z = boot_log2_mean / boot_log2_se
  boot_log2_pValue = 2*stats::pnorm(q = abs(boot_log2_z), lower.tail = FALSE)

  # boot_mean <- 2^boot_log2_mean
  # boot_se <- apply(boot[], 1, function(x) stats::sd(2^x, na.rm = TRUE))
  # boot_z = boot_mean / boot_se
  # boot_pValue = 2*stats::pnorm(q = abs(boot_z), lower.tail = FALSE)

  #create results dataframe
  df_results <- data.frame(non_boot_means,
                           boot_log2_se = boot_log2_se,
                           boot_log2_z = boot_log2_z,
                           boot_log2_pValue = boot_log2_pValue,
                           valid_boot_samples = N_valid_boot_resamples)

  #print results
  print(paste0("On average, ", round(mean(N_valid_boot_resamples),2), " out of ",  n_boot, " bootstrap resamples produced valid mean diagnosticity scores"))
  print(df_results)

  #if no failure up to this point, update the results list with computed values
  results$df_results <- df_results
  results$boot <- boot

  return(results)
}
