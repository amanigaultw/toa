#' decouples multi-gene rows
#'
#' produces a data frame containing one row per gene symbols, where expression values
#' are identical for previously coupled symbols.
#'
#' @param df a data frame or matrix containing gene symbols and gene expression values.
#' @param sep a string indicating how multi-gene symbols are separated (e.g., "///").
#' @param gene_col a numeric value indicating the position of the gene symbol column
#' within \code{df}; the first column is used by default.
#' @param header bool indicating whether the first row of data is a header; TRUE by default.
#' @param progress bool indicating whether to display function progress.
#' @return a data frame containing one row per gene symbols, where expression values are
#' identical for previously coupled symbols.
#' @examples
#' \dontrun{
#' #load example data
#' data("HumanCD14CD16NegVsPosVsDC_3Reps_Ziegler")
#' #decouple_genes
#' decoupled_data <- decouple_genes(HumanCD14CD16NegVsPosVsDC_3Reps_Ziegler)
#' }
#' @export
decouple_genes <- function(df, sep = "///", gene_col = 1, header = TRUE, progress = TRUE){

  df <- as.data.frame(df)

  #store header row
  if(header){
    temp_row <- df[1,]
    df <- df[-1,]
  }

  #set expression columns to numeric
  df[,-gene_col] <- sapply(df[,-gene_col], as.numeric)
  #get rid of empty gene symbol rows
  df <- subset(df, df[,gene_col] != "")

  #decouple multi-gene rows
  n <- length(unlist(strsplit(df[,gene_col], sep, fixed = TRUE)))
  temp <- vector("list", n)
  for(i in 1:nrow(df)){
    symbol <- df[i,gene_col]

    if(grepl(sep, symbol)){
      symbols <- unlist(strsplit(symbol, sep, fixed = TRUE)) #generate vector of gene symbols
      symbols <- gsub(" ", "", symbols) #get rid of spaces in gene symbol names
      m <- matrix(nrow = length(symbols), ncol = ncol(df))
      for(j in 1:length(symbols)){
        m[j,] <- c(symbols[j], t(df[i,-gene_col]))
      }
      temp[[i]] <- m
    }else{
      temp[[i]] <- as.matrix(df[i,])
    }

    if(progress == TRUE){
      p <- i / nrow(df) * 100
      svMisc::progress(p)
    }
  }

  #generate a dataframe from the list
  df_temp <- do.call(rbind, temp)

  #clean gene symbols
  df_temp[,1] <- trimws(toupper(df_temp[,1]))

  #restore header
  if(header){
    df_temp <- rbind(temp_row, df_temp)
    row.names(df_temp) <- NULL
  }

  return(df_temp)
}


