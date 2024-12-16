#' Check if two data frames have matching sample names and order
#'
#' This function checks if two data frames have matching rownames and row order. If the rownames
#' match but are in a different order, the function will reorder both data frames to match.
#'
#' @param df1 A data frame with data, where rows represent samples
#' @param df2 A data frame with data, where rows represent samples
#' @return A list with both df1 and df2 data in the same order and with matching rownames
#' @examples
#' data1 <- data.frame(A=c(1,2,3), B=c(4,5,6), row.names=c("S1", "S2", "S3"))
#' data2 <- data.frame(X=c(7,8,9), Y=c(10,11,12), row.names=c("S2", "S1", "S3"))
#' df_list <- CheckSampleOrder(data1, data2)
#' data1_reordered <- df_list[[1]]
#' data2_reordered <- df_list[[2]]
#' @export
CheckSampleOrder <- function(df1, df2){
  same_order <- all(rownames(df1) == rownames(df2))
  if(!same_order){
    cat("Sample names do not match. Samples will be sorted to match.")
    df1 <- df1[order(rownames(df1)),]
    df2 <- df2[order(rownames(df2)),]
  } else{
    cat("Samples match")
  }
  return(list(df1, df2))
}


