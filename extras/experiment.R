library(toa)

data("Chang")
data("epith_mesen_ref")



scores <- toa::compute_diagnosticity_scores(Dataframe = epith_mesen_ref,
                                            gene_col = 1,
                                            Cell_of_interest_col = 2:11,
                                            Reference_cells_col = 12:21,
                                            logZeroSubstituteExpressionValue = .001,
                                            na.rm = TRUE,
                                            gene_to_upper = TRUE)

test1 <- toa::gene_TOA_function(Analysis_Dataframe = Chang,
                                predictor_col = "stress",
                                TOA_Reference_Dataframe = scores,
                                covariate_cols = NULL,
                                gene_cols = NULL,
                                ref_gene_col = "gene",
                                foldThreshDEG = 1.25)
test1$df.means

test2 <- toa::bootstrap_gene_TOA_function(Bootstrap_Samples = 200,
                                          Analysis_Dataframe = Chang,
                                          predictor_col = "stress",
                                          TOA_Reference_Dataframe = scores,
                                          covariate_cols = NULL,
                                          gene_cols = NULL,
                                          ref_gene_col = "gene",
                                          foldThreshDEG = 1.25)

test2$results


