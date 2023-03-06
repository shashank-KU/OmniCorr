#' Calculate correlation between two data frames
#'
#' @description
#' This function calculates the correlation matrix between the columns of two data frames.
#'
#' @param df1 A data frame containing the first set of variables to be correlated. The samples should be in the rows and the MEs in the columns.
#' @param df2 A data frame containing the second set of variables to be correlated. The samples should be in the rows and the MEs in the columns.
#' @param method (optional) The correlation method to be used. Default is "pearson".
#' @param adjust_method (optional) The method for adjusting the p-values for multiple testing. Default is "fdr".
#'
#' @return A list containing the correlation results
#' @return$correlation Matrix of correlation coefficients.
#' @return$p_value Matrix of p-values for the correlation coefficients.
#' @return$p_value_adj Matrix of adjusted p-values for the correlation coefficients.
#' @return$signif_matrix Matrix indicating the significance level of the correlations.
#'
#' @examples
#' # Calculate correlation between two data frames using default parameters
#' result_list <- calculate_correlations(df1 = Transcriptomics, df2 = Metagenomics)
#'
#' @export
calculate_correlations <- function(df1, df2, method="pearson", adjust_method="fdr") {
  # Calculate correlation matrix
  cor_mat <- cor(df1, df2, method=method)

  # Calculate p-values and adjust for multiple testing
  p_val_mat <- matrix(ncol=ncol(df2), nrow=ncol(df1),
                      dimnames=list(colnames(df1), colnames(df2)))
  corr_val_mat <- p_val_mat

  for (i in 1:ncol(df1)) {
    for (j in 1:ncol(df2)) {
      cor_res <- cor.test(df1[,i], df2[,j], method=method)
      p_val_mat[i,j] <- cor_res$p.value
      corr_val_mat[i,j] <- cor_res$estimate
    }
  }

  p_val_adj_mat <- p.adjust(p_val_mat, method=adjust_method)

  # Add significance level.
  # One star means a p-value of less than 0.05;
  # Two stars is less than 0.01, and three, is less than 0.001.
  signif_matrix <- rep("", length(p_val_mat))
  three_star <- which(p_val_mat <= 0.001)
  signif_matrix[three_star] <- "***"
  two_star <- which((p_val_mat <= 0.01) & (p_val_mat > 0.001))
  signif_matrix[two_star] <- "**"
  one_star <- which((p_val_mat <= 0.05) & (p_val_mat > 0.01))
  signif_matrix[one_star] <- "*"
  dim(signif_matrix) <- dim(p_val_mat) # Give textMatrix the correct dimensions

  # Collect all results into a list
  results <- list(p_value=p_val_mat,
                  p_value_adj=p_val_adj_mat,
                  signif_matrix=signif_matrix,
                  correlation=corr_val_mat)

  return(results)
}




