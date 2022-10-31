test_that("toa produces the same output as Steve's Java program on TAU trial data", {

  #load example data
  data("TAU_Trials3_Gene_CPM_Log2")
  data("TAUTrials2022BC_Intervention_Rm1BadQC_RmB14")
  data("HumanM1M2_3Reps_Martinez")

  popGenes <- read.table(test_path("fixtures", "BootstrapMicroarrayAnalysis - TAU2022BC - Intervention RmB14 - EMT Transfac 2.00 - popGenes.txt"))
  popGenes <- as.character(popGenes)
  popGenes[which(popGenes == "TRUE")] <- "T" #deal with annoying re-coding of gene name by Steve's program.

  #get DEG
  DEG_result <- get_DEG(expression_data = TAU_Trials3_Gene_CPM_Log2,
                        regressor_matrix = TAUTrials2022BC_Intervention_Rm1BadQC_RmB14,
                        foldThreshDEG = 2,
                        verbose = F)
  #get toa
  toa_result <- toa(DEG_result = DEG_result,
                    ref = HumanM1M2_3Reps_Martinez,
                    type_1_cols = 2:4,
                    type_2_cols = 5:7,
                    verbose = F)

  expect_equal(round(toa_result$df_results, 4)[1,1], .0844)
  expect_equal(round(toa_result$df_results, 4)[2,1], .3036)
  expect_equal(round(toa_result$df_results, 4)[3,1], .1666)
  expect_equal(round(toa_result$df_results, 4)[4,1], .0841)
  expect_equal(round(toa_result$df_results, 4)[5,1], .3596)
  expect_equal(round(toa_result$df_results, 4)[6,1], .1874)

  expect_equal(sort(toa_result$pop_gene_symbols), sort(popGenes))
})
