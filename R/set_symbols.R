#' reformat gene symbols
#'
#' reformats gene symbols such that space characters are omitted and all remaining
#' characters are upper case.
#'
#' @param symbols a character vector of gene symbols.
#' @return a reformatted character vector of gene symbols.
#' @examples
#' #load example data
#' data("epith_mesen_ref_raw")
#' #clean symbols
#' epith_mesen_ref_raw$gene <- set_symbols(epith_mesen_ref_raw$gene)
#' @export
set_symbols <- function(symbols){
  return(trimws(toupper(symbols)))
}
