#' Calculate correlation between two data frames
#'
#' @description
#' This function calculates the correlation matrix between the columns of two data frames.
#'
#' @param df1 A data frame or matrix containing the first set of variables to be correlated. Rows are samples, columns are variables (e.g., Genes/Proteins/Metabolites/ASVs).
#' @param df2 A data frame or matrix containing the second set of variables to be correlated. Rows are samples, columns are variables (e.g., Genes/Proteins/Metabolites/ASVs).
#' @param method (optional) The correlation method to be used. Default is "pearson". Other options include "spearman" and "kendall".
#' @param adjust_method (optional) The method for adjusting the p-values for multiple testing. Default is "fdr". If "none", no adjustment is applied. Supports all methods from \code{p.adjust.methods}.
#' @param use (optional) A character string specifying the handling of missing data. Default is "all.obs". Currently supports "all.obs" and "pairwise.complete.obs". Abbreviations are not allowed in this version for clarity.
#' @param show_significance (optional) A character string indicating how to represent significance in the signif_matrix. Default is "stars". Possible values are "stars" (significance stars), "p_value" (adjusted p-values), or "correlation" (correlation coefficients).
#'
#' @return A list containing the correlation results:
#' \item{correlation}{Matrix of correlation coefficients (rounded to 3 decimal places).}
#' \item{p_value}{Matrix of raw p-values for the correlation coefficients.}
#' \item{p_value_adj}{Matrix of adjusted p-values for the correlation coefficients (or raw if no adjustment).}
#' \item{signif_matrix}{Matrix indicating significance: stars ("***" for p_adj <= 0.001, "**" for <= 0.01, "*" for <= 0.05), adjusted p-values, or correlation coefficients, based on \code{show_significance}.}
#'
#' @details
#' For "pearson" and "spearman", efficient computation is used via \code{WGCNA::corAndPvalue}. For "kendall", falls back to a loop for p-values, which may be slow for large datasets.
#' Assumes columns are numeric; no explicit checks are performed for performance, but non-numeric data will cause errors.
#' For \code{use = "all.obs"}, an error is thrown if any NAs are present in the combined data.
#' For \code{use = "pairwise.complete.obs"}, NAs are handled per pair.
#' Other \code{use} options from base \code{cor} are not supported in this version.
#'
#' @examples
#' # Calculate correlation between two data frames using default parameters and show significance as p-values
#' result_list <- calculate_correlations(df1 = Transcriptomics, df2 = Metagenomics, show_significance = "p_value")
#'
#' @export
calculate_correlations <- function(df1, df2, method = "pearson", adjust_method = "fdr", use = "all.obs", show_significance = "stars") {
  # Load required package
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' must be installed to use this function.")
  }

  # Validate inputs
  if (!is.data.frame(df1) && !is.matrix(df1)) stop("df1 must be a data frame or matrix.")
  if (!is.data.frame(df2) && !is.matrix(df2)) stop("df2 must be a data frame or matrix.")
  if (nrow(df1) != nrow(df2)) stop("df1 and df2 must have the same number of rows (samples).")
  if (ncol(df1) == 0 || ncol(df2) == 0) stop("df1 and df2 must have at least one column.")

  valid_methods <- c("pearson", "spearman", "kendall")
  method <- match.arg(method, valid_methods)

  valid_uses <- c("all.obs", "pairwise.complete.obs")
  use <- match.arg(use, valid_uses)

  valid_adjust <- c("none", p.adjust.methods)
  adjust_method <- match.arg(adjust_method, valid_adjust)

  valid_show <- c("stars", "p_value", "correlation")
  show_significance <- match.arg(show_significance, valid_show)

  # Handle NA behavior for consistency
  combined_data <- cbind(df1, df2)
  if (use == "all.obs" && !all(complete.cases(combined_data))) {
    stop("Data contains NAs, but use = 'all.obs' requires complete observations.")
  }

  # Compute correlation and p-values
  if (method %in% c("pearson", "spearman")) {
    res <- WGCNA::corAndPvalue(df1, df2, method = method, use = use)
    cor_mat <- res$cor
    p_val_mat <- res$p
  } else {  # kendall
    cor_mat <- cor(df1, df2, method = method, use = use)
    p_val_mat <- matrix(NA_real_, nrow = ncol(df1), ncol = ncol(df2),
                        dimnames = list(colnames(df1), colnames(df2)))
    for (i in 1:ncol(df1)) {
      for (j in 1:ncol(df2)) {
        x <- df1[, i]
        y <- df2[, j]
        # cor.test handles NAs pairwise; since use="all.obs" already checked globally, it's fine
        p_val_mat[i, j] <- cor.test(x, y, method = method)$p.value
      }
    }
  }

  # Adjust p-values
  if (adjust_method != "none") {
    p_val_adj_mat <- matrix(p.adjust(p_val_mat, method = adjust_method),
                            nrow = nrow(p_val_mat), ncol = ncol(p_val_mat),
                            dimnames = dimnames(p_val_mat))
  } else {
    p_val_adj_mat <- p_val_mat
  }

  # Create signif_matrix
  signif_matrix <- matrix(NA_real_, nrow = nrow(p_val_mat), ncol = ncol(p_val_mat),
                          dimnames = dimnames(p_val_mat))
  if (show_significance == "stars") {
    signif_matrix_char <- matrix("", nrow = nrow(p_val_mat), ncol = ncol(p_val_mat),
                                 dimnames = dimnames(p_val_mat))
    signif_matrix_char[p_val_adj_mat <= 0.001] <- "***"
    signif_matrix_char[p_val_adj_mat <= 0.01 & p_val_adj_mat > 0.001] <- "**"
    signif_matrix_char[p_val_adj_mat <= 0.05 & p_val_adj_mat > 0.01] <- "*"
    signif_matrix <- signif_matrix_char  # Character matrix for stars
  } else if (show_significance == "p_value") {
    signif_matrix[] <- p_val_adj_mat
  } else if (show_significance == "correlation") {
    signif_matrix[] <- cor_mat
  }

  # Round correlation matrix for return (to match original behavior)
  correlation_rounded <- round(cor_mat, 3)

  # Return results
  results <- list(
    correlation = correlation_rounded,
    p_value = p_val_mat,
    p_value_adj = p_val_adj_mat,
    signif_matrix = signif_matrix
  )
  return(results)
}
