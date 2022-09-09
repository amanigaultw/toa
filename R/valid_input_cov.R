valid_input_cov <- function(cov){

  if(!is.null(cov)){
    sds <- apply(as.matrix(cov), 2, function(x) sd(x, na.rm = TRUE))

    if(sum(is.na(cov)) > 0 | any(sds == 0)){
      print("invalid cov input")
      return(FALSE)
    }
  }

  return(TRUE)
}
