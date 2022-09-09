get_means <- function(df.DEG.ref.matched, df.DEG.ref.matched.expressed){

  #compute population diagnosticity score
  cellDiagnosticityPopulationMeanScoresLinear <- mean(df.DEG.ref.matched$DiagnosticityScoresLinear)
  cellDiagnosticityPopulationMeanScoresLog <- mean(df.DEG.ref.matched$DiagnosticityScoresLog)

  # aggregate by DEG (i.e., flagged as up vs. down regulated)
  df.temp.means <- aggregate(df.DEG.ref.matched.expressed[, 4:5], list(df.DEG.ref.matched.expressed$DEG), mean)

  #compute results
  up_reg_linear = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
  down_reg_linear = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLinear"] - cellDiagnosticityPopulationMeanScoresLinear
  reg_linear = mean(df.DEG.ref.matched.expressed$DiagnosticityScoresLinear) - cellDiagnosticityPopulationMeanScoresLinear
  up_reg_log = df.temp.means[df.temp.means[,1] == "1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
  down_reg_log = df.temp.means[df.temp.means[,1] == "-1", "DiagnosticityScoresLog"] - cellDiagnosticityPopulationMeanScoresLog
  reg_log = mean(df.DEG.ref.matched.expressed$DiagnosticityScoresLog) - cellDiagnosticityPopulationMeanScoresLog

  return(list(up_reg_linear = up_reg_linear,
              down_reg_linear = down_reg_linear,
              reg_linear = reg_linear,
              up_reg_log = up_reg_log,
              down_reg_log = down_reg_log,
              reg_log = reg_log))
}
