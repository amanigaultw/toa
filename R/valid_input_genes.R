valid_input_genes <- function(genes){

  if(sum(is.na(genes)) > 0){
    print("invalid gene input")
    return(FALSE)
  }

  return(TRUE)
}
