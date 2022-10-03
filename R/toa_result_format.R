#' Format TOA results
#'
#' formats results from the bootstrapped transcript origin analysis
#' function to ease comparison with past/other work.
#'
#' @param toa_boot_result a toa_boot result object produced using \code{toa_boot()}.#' @return toa returns an object of class "toa"
#' @return a formated results matrix.
#' @examples
#' \dontrun{
#' #load example data
#' data("TAU_Trials3_Gene_CPM_Log2")
#' data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
#' data("HumanM1M2_3Reps_Martinez")
#'
#' #get DEG
#' DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
#'                       regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
#'                       foldThreshDEG = 2)
#'
#' #get toa
#' toa_result <- toa(DEG_result = DEG_result,
#'                   ref = HumanM1M2_3Reps_Martinez,
#'                   type_1_cols = 2:4,
#'                   type_2_cols = 5:7)
#'
#' #get bootstrapped stats
#' toa_boot_result <- toa_boot(toa_result)
#'
#' #generate a formated toa results table
#' pretty_results <- toa_result_format(toa_boot_result)
#' View(pretty_results)
#' }
#' @export
toa_result_format <- function(toa_boot_result){

steve_labels_1 <- c(names(toa_boot_result$inputs$toa$inputs$DEG$inputs$regressor_matrix)[2],
                    "",
                    "Linear diagnosticity scores",
                    "Cell type",
                    "Mean diagnosticity score",
                    "Mean Bootstrap",
                    "SE Bootstrap",
                    "p SE Bootstrap",
                    "p quantile Bootstrap",
                    paste("n genes | fold-dif > ", format(toa_boot_result$inputs$toa$inputs$DEG$inputs$foldThreshDEG, nsmall = 2)),
                    "Log diagnosticity scores",
                    "Cell type",
                    "Mean diagnosticity score",
                    "Mean Bootstrap", "SE Bootstrap",
                    "p SE Bootstrap",
                    "p quantile Bootstrap",
                    paste("n genes | fold-dif > ", format(toa_boot_result$inputs$toa$inputs$DEG$inputs$foldThreshDEG, nsmall = 2)))

steve_labels_2 <- c("Linear diagnosticity scores",
                    "Upregulated",
                    "",
                    "",
                    "Downregulated",
                    "",
                    "",
                    "Regulated",
                    "",
                    "",
                    "")

steve_labels_3 <- c("Cell type",
                    "cell type 1",
                    "cell type 2",
                    "",
                    "cell type 1",
                    "cell type 2",
                    "",
                    "cell type 1",
                    "cell type 2",
                    "",
                    "")


table <- matrix("", nrow = 18, ncol = 11)
table[3, ] <- steve_labels_2
table[4, ] <- steve_labels_3
table[11, ] <- steve_labels_2
table[12, ] <- steve_labels_3
table[,1] <- steve_labels_1

#declare some variable before using list2env()
mean_toa_upreg_linear <- mean_toa_downreg_linear <- mean_toa_reg_linear <- NULL
mean_toa_upreg_log <- mean_toa_downreg_log <- mean_toa_reg_log <- NULL
boot_mean_toa_upreg_linear <- boot_mean_toa_downreg_linear <- boot_mean_toa_reg_linear <- NULL
boot_mean_toa_upreg_log <- boot_mean_toa_downreg_log <- boot_mean_toa_reg_log <- NULL
boot_se_toa_upreg_linear <- boot_se_toa_downreg_linear <- boot_se_toa_reg_linear <- NULL
boot_se_toa_upreg_log <- boot_se_toa_downreg_log <- boot_se_toa_reg_log <- NULL
p_value_upreg_linear <- p_value_downreg_linear <- p_value_reg_linear <- NULL
p_value_upreg_log <- p_value_downreg_log <- p_value_reg_log <- NULL

list2env(toa_boot_result$inputs$toa$toa_means, envir=environment())
list2env(toa_boot_result$boot_stats$boot_mean, envir=environment())
list2env(toa_boot_result$boot_stats$boot_se, envir=environment())
list2env(toa_boot_result$boot_stats$boot_p_value, envir=environment())
df_DEG <- toa_boot_result$inputs$toa$inputs$DEG$df_DEG

table[5,3] <- round(mean_toa_upreg_linear, 4)
table[5,6] <- round(mean_toa_downreg_linear, 4)
table[5,9] <- round(mean_toa_reg_linear, 4)
table[13,3] <- round(mean_toa_upreg_log, 4)
table[13,6] <- round(mean_toa_downreg_log, 4)
table[13,9] <- round(mean_toa_reg_log, 4)

table[6,3] <- round(boot_mean_toa_upreg_linear, 4)
table[6,6] <- round(boot_mean_toa_downreg_linear, 4)
table[6,9] <- round(boot_mean_toa_reg_linear, 4)
table[14,3] <- round(boot_mean_toa_upreg_log, 4)
table[14,6] <- round(boot_mean_toa_downreg_log, 4)
table[14,9] <- round(boot_mean_toa_reg_log, 4)

table[7,3] <- round(boot_se_toa_upreg_linear, 4)
table[7,6] <- round(boot_se_toa_downreg_linear, 4)
table[7,9] <- round(boot_se_toa_reg_linear, 4)
table[15,3] <- round(boot_se_toa_upreg_log, 4)
table[15,6] <- round(boot_se_toa_downreg_log, 4)
table[15,9] <- round(boot_se_toa_reg_log, 4)

table[8,3] <- round(p_value_upreg_linear, 4)
table[8,6] <- round(p_value_downreg_linear, 4)
table[8,9] <- round(p_value_reg_linear, 4)
table[16,3] <- round(p_value_upreg_log, 4)
table[16,6] <- round(p_value_downreg_log, 4)
table[16,9] <- round(p_value_reg_log, 4)

table[10,3] <- table(df_DEG$DEG)[names(table(df_DEG$DEG)) == 1]
table[10,6] <- table(df_DEG$DEG)[names(table(df_DEG$DEG)) == -1]
table[10,9] <- sum(table(df_DEG$DEG)[names(table(df_DEG$DEG)) != 0])
table[18,3] <- table(df_DEG$DEG)[names(table(df_DEG$DEG)) == 1]
table[18,6] <- table(df_DEG$DEG)[names(table(df_DEG$DEG)) == -1]
table[18,9] <- sum(table(df_DEG$DEG)[names(table(df_DEG$DEG)) != 0])

table[5,2] <- -as.numeric(table[5,3])
table[5,5] <- -as.numeric(table[5,6])
table[5,8] <- -as.numeric(table[5,9])
table[13,2] <- -as.numeric(table[13,3])
table[13,5] <- -as.numeric(table[13,6])
table[13,8] <- -as.numeric(table[13,9])

table[6,2] <- -as.numeric(table[6,3])
table[6,5] <- -as.numeric(table[6,6])
table[6,8] <- -as.numeric(table[6,9])
table[14,2] <- -as.numeric(table[14,3])
table[14,5] <- -as.numeric(table[14,6])
table[14,8] <- -as.numeric(table[14,9])

table[7,2] <- as.numeric(table[7,3])
table[7,5] <- as.numeric(table[7,6])
table[7,8] <- as.numeric(table[7,9])
table[15,2] <- as.numeric(table[15,3])
table[15,5] <- as.numeric(table[15,6])
table[15,8] <- as.numeric(table[15,9])

table[8,2] <- 1 - as.numeric(table[8,3])
table[8,5] <- 1 - as.numeric(table[8,6])
table[8,8] <- 1 - as.numeric(table[8,9])
table[16,2] <- 1 - as.numeric(table[16,3])
table[16,5] <- 1 - as.numeric(table[16,6])
table[16,8] <- 1 - as.numeric(table[16,9])

table[10,2] <- table[10,3]
table[10,5] <- table[10,6]
table[10,8] <- table[10,9]
table[18,2] <- table[18,3]
table[18,5] <- table[18,6]
table[18,8] <- table[18,9]

return(table)
}
