valid_input_x <- function(x){

  if(sum(is.na(x)) > 0 | sd(x) == 0) {
    print("invalid x input")
    return(FALSE)
  }

  return(TRUE)
}
