valid_input_genes <- function(genes){

  if(sum(is.na(genes)) > 0){
    warning("invalid gene input")
    return(FALSE)
  }

  # sds <- apply(as.matrix(genes), 2, function(x) sd(x, na.rm = TRUE))
  # if(any(sds) == 0){
  #   print("invalid gene input")
  #   return(FALSE)
  # }

  return(TRUE)
}
