test_that("get_DEG works", {

  data("Chang")
  DEG_result <- get_DEG(x = Chang$stress,
                        genes = subset(Chang, select = -stress),
                        foldThreshDEG = 1.25)
  temp <- table(DEG_result$DEG)
  names(temp) <- NULL
  downreg_count <- temp[1]
  upreg_count <- temp[3]

  expect_equal(downreg_count, 123)
  expect_equal(upreg_count, 227)
})
