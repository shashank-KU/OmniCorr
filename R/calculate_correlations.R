#' Calculate correlation between two data frames
#'
#' @description
#' This function calculates the correlation matrix between the columns of two data frames.
#'
#' @param df1 A data frame or matrix containing the first set of variables to be correlated.
#'   Rows are samples, columns are variables (e.g., genes, proteins, metabolites).
#' @param df2 A data frame or matrix containing the second set of variables to be correlated.
#'   Rows are samples, columns are variables.
#' @param method The correlation method to be used.
#'   Default is "pearson". Other options include "spearman" and "kendall".
#' @param adjust_method Method for adjusting the p-values for multiple testing.
#'   Default is "fdr". If "none", no adjustment is applied.
#'   Supports all methods from \code{p.adjust.methods}.
#' @param use Handling of missing data.
#'   One of "all.obs" or "pairwise.complete.obs".
#' @param show_significance How to represent significance in the output.
#'   One of "stars", "p_value", or "correlation".
#'
#' @return A list containing:
#' \item{correlation}{Matrix of correlation coefficients (rounded to 3 decimals).}
#' \item{p_value}{Matrix of raw p-values.}
#' \item{p_value_adj}{Matrix of adjusted p-values.}
#' \item{signif_matrix}{Matrix of significance stars, p-values, or correlations.}
#'
#' @details
#' Pearson and Spearman correlations are computed using
#' \code{WGCNA::corAndPvalue} for efficiency.
#' Kendall correlation uses \code{cor.test} in a loop and may be slow for large datasets.
#'
#' @examples
#' data(Transcriptomics)
#' data(Metagenomics)
#'
#' result_list <- calculate_correlations(
#'   df1 = Transcriptomics,
#'   df2 = Metagenomics,
#'   show_significance = "p_value"
#' )
#'
#' @importFrom stats cor cor.test complete.cases p.adjust p.adjust.methods
#' @export
calculate_correlations <- function(
    df1,
    df2,
    method = "pearson",
    adjust_method = "fdr",
    use = "all.obs",
    show_significance = "stars"
) {
  
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' must be installed.")
  }
  
  if (!is.data.frame(df1) && !is.matrix(df1))
    stop("df1 must be a data frame or matrix.")
  if (!is.data.frame(df2) && !is.matrix(df2))
    stop("df2 must be a data frame or matrix.")
  if (nrow(df1) != nrow(df2))
    stop("df1 and df2 must have the same number of rows (samples).")
  
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  use <- match.arg(use, c("all.obs", "pairwise.complete.obs"))
  adjust_method <- match.arg(adjust_method, c("none", stats::p.adjust.methods))
  show_significance <- match.arg(
    show_significance,
    c("stars", "p_value", "correlation")
  )
  
  combined_data <- cbind(df1, df2)
  if (use == "all.obs" && !all(stats::complete.cases(combined_data))) {
    stop("Missing values present but use = 'all.obs'.")
  }
  
  if (method %in% c("pearson", "spearman")) {
    res <- WGCNA::corAndPvalue(df1, df2, method = method, use = use)
    cor_mat <- res$cor
    p_val_mat <- res$p
  } else {
    cor_mat <- stats::cor(df1, df2, method = method, use = use)
    p_val_mat <- matrix(
      NA_real_,
      nrow = ncol(df1),
      ncol = ncol(df2),
      dimnames = list(colnames(df1), colnames(df2))
    )
    for (i in seq_len(ncol(df1))) {
      for (j in seq_len(ncol(df2))) {
        p_val_mat[i, j] <- stats::cor.test(
          df1[, i], df2[, j],
          method = method
        )$p.value
      }
    }
  }
  
  if (adjust_method != "none") {
    p_val_adj_mat <- matrix(
      stats::p.adjust(p_val_mat, method = adjust_method),
      nrow = nrow(p_val_mat),
      dimnames = dimnames(p_val_mat)
    )
  } else {
    p_val_adj_mat <- p_val_mat
  }
  
  if (show_significance == "stars") {
    signif_matrix <- matrix(
      "",
      nrow = nrow(p_val_adj_mat),
      ncol = ncol(p_val_adj_mat),
      dimnames = dimnames(p_val_adj_mat)
    )
    signif_matrix[p_val_adj_mat <= 0.001] <- "***"
    signif_matrix[p_val_adj_mat <= 0.01 & p_val_adj_mat > 0.001] <- "**"
    signif_matrix[p_val_adj_mat <= 0.05 & p_val_adj_mat > 0.01] <- "*"
  } else if (show_significance == "p_value") {
    signif_matrix <- p_val_adj_mat
  } else {
    signif_matrix <- cor_mat
  }
  
  list(
    correlation = round(cor_mat, 3),
    p_value = p_val_mat,
    p_value_adj = p_val_adj_mat,
    signif_matrix = signif_matrix
  )
}
