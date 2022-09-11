.valid_input_toa_ref <- function(toa_ref){

  if(!all.equal(colnames(toa_ref), c("gene", "DiagnosticityScoresLinear", "DiagnosticityScoresLog")) | !nrow(toa_ref) > 0) {
    warning("invalid toa_ref input")
    return(FALSE)
  }

  return(TRUE)
}
