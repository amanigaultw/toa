test_that("get_DEG works", {

  #load example data
  data("TAU_Trials3_Gene_CPM_Log2")
  data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")

  #get DEG
  DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
                        regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
                        foldThreshDEG = 2)

  temp <- table(DEG_result$df_DEG$DEG)
  names(temp) <- NULL
  downreg_count <- temp[1]
  upreg_count <- temp[3]

  expect_equal(downreg_count, 1208)
  expect_equal(upreg_count, 1962)
})
