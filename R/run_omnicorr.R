#' Run OmniCorr Multi-Omics Integration (Production Version)
#'
#' @description
#' Integrates multiple omics layers by computing correlations between a
#' reference layer and comparison layers, optionally including metadata.
#' Produces a multi-panel integrated heatmap using ComplexHeatmap.
#'
#' @inheritParams calculate_correlations
#'
#' @param reference_layer Numeric matrix/data.frame (samples x features).
#' @param comparison_layers Named list of numeric matrices/data.frames.
#' @param metadata Optional numeric matrix/data.frame.
#'
#' @param cluster_reference_rows Logical.
#' @param cluster_reference_columns Logical.
#' @param cluster_comparison_rows Logical.
#' @param cluster_comparison_columns Logical.
#' @param cluster_metadata_rows Logical.
#' @param cluster_metadata_columns Logical.
#'
#' @param show_row_names Logical.
#' @param show_column_names Logical.
#' @param font_size Numeric.
#' @param star_font_size Numeric.
#'
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly returns a ComplexHeatmap::HeatmapList object.
#'
#' @export
run_omnicorr <- function(
    reference_layer,
    comparison_layers,
    metadata = NULL,
    method = "pearson",
    adjust_method = "fdr",
    use = "all.obs",
    show_significance = "stars",
    cluster_reference_rows = TRUE,
    cluster_reference_columns = FALSE,
    cluster_comparison_rows = FALSE,
    cluster_comparison_columns = FALSE,
    cluster_metadata_rows = FALSE,
    cluster_metadata_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    font_size = 8,
    star_font_size = 8,
    verbose = TRUE
) {
  
  # ----------------------------
  # Required packages
  # ----------------------------
  required_pkgs <- c("WGCNA", "ComplexHeatmap", "circlize", "RColorBrewer", "grid")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' must be installed.", pkg))
    }
  }
  
  # ----------------------------
  # Argument validation
  # ----------------------------
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  use <- match.arg(use, c("all.obs", "pairwise.complete.obs"))
  adjust_method <- match.arg(adjust_method, c("none", stats::p.adjust.methods))
  show_significance <- match.arg(show_significance,
                                 c("stars", "p_value", "correlation"))
  
  if (!is.matrix(reference_layer) && !is.data.frame(reference_layer))
    stop("reference_layer must be a matrix or data.frame.")
  
  reference_layer <- as.matrix(reference_layer)
  
  if (!all(apply(reference_layer, 2, is.numeric)))
    stop("reference_layer must contain only numeric values.")
  
  if (is.null(colnames(reference_layer)))
    stop("reference_layer must have column names.")
  
  if (is.null(rownames(reference_layer)))
    stop("reference_layer must have row names (sample IDs).")
  
  if (!is.list(comparison_layers))
    stop("comparison_layers must be a named list.")
  
  if (is.null(names(comparison_layers)) ||
      any(names(comparison_layers) == ""))
    stop("comparison_layers must be a named list.")
  
  # Validate comparison layers
  for (nm in names(comparison_layers)) {
    
    layer <- as.matrix(comparison_layers[[nm]])
    
    if (!all(apply(layer, 2, is.numeric)))
      stop(sprintf("Layer '%s' contains non-numeric values.", nm))
    
    if (!identical(rownames(reference_layer), rownames(layer)))
      stop(sprintf("Sample order mismatch in layer '%s'.", nm))
    
    comparison_layers[[nm]] <- layer
  }
  
  # Validate metadata
  if (!is.null(metadata)) {
    
    metadata <- as.matrix(metadata)
    
    if (!all(apply(metadata, 2, is.numeric)))
      stop("metadata must contain only numeric values.")
    
    if (!identical(rownames(reference_layer), rownames(metadata)))
      stop("Sample order mismatch between reference_layer and metadata.")
  }
  
  # ----------------------------
  # Progress header
  # ----------------------------
  if (verbose) {
    message("==================================================")
    message("🔷 OmniCorr Multi-Omics Integration")
    message("==================================================")
    message("Samples: ", nrow(reference_layer))
    message("Reference Features: ", ncol(reference_layer))
    message("✔ Sample alignment verified")
  }
  
  # ----------------------------
  # Reference clustering
  # ----------------------------
  if (cluster_reference_rows) {
    
    if (verbose) message("Clustering reference rows...")
    
    cor_ref <- WGCNA::bicor(reference_layer, maxPOutliers = 0.05)
    cor_ref[is.na(cor_ref)] <- 0
    
    dendro <- stats::hclust(stats::as.dist(1 - cor_ref))
    ordered_features <- colnames(reference_layer)[dendro$order]
    
  } else {
    
    ordered_features <- colnames(reference_layer)
    
    if (verbose)
      message("Reference row clustering skipped")
  }
  
  # ----------------------------
  # Color functions
  # ----------------------------
  rdBu_colors <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(6, "RdBu"))
  )(51)
  
  cor_col_fun <- circlize::colorRamp2(
    seq(-1, 1, length.out = 51),
    rdBu_colors
  )
  
  expr_range <- range(reference_layer, na.rm = TRUE)
  
  expr_col_fun <- circlize::colorRamp2(
    seq(expr_range[1], expr_range[2], length.out = 51),
    rdBu_colors
  )
  
  # ----------------------------
  # Heatmap helper
  # ----------------------------
  create_cor_ht <- function(mat, signif_mat, title_text,
                            cluster_rows_flag, cluster_cols_flag) {
    
    safe_width <- min(ncol(mat) * 3, 250)
    
    ComplexHeatmap::Heatmap(
      mat,
      name = "Correlation",
      col = cor_col_fun,
      cluster_rows = cluster_rows_flag,
      cluster_columns = cluster_cols_flag,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = grid::gpar(fontsize = font_size),
      column_names_gp = grid::gpar(fontsize = font_size),
      column_title = title_text,
      width = grid::unit(safe_width, "mm"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        
        if (show_significance == "stars") {
          grid::grid.text(signif_mat[i, j], x, y,
                          gp = grid::gpar(fontsize = star_font_size))
        } else if (show_significance == "p_value") {
          grid::grid.text(sprintf("%.2g", signif_mat[i, j]), x, y,
                          gp = grid::gpar(fontsize = star_font_size))
        } else {
          grid::grid.text(sprintf("%.2f", signif_mat[i, j]), x, y,
                          gp = grid::gpar(fontsize = star_font_size))
        }
      }
    )
  }
  
  ht_list <- NULL
  
  # ----------------------------
  # Metadata
  # ----------------------------
  if (!is.null(metadata)) {
    
    if (verbose) message("Computing metadata correlations...")
    
    meta_res <- calculate_correlations(
      reference_layer,
      metadata,
      method,
      adjust_method,
      "pairwise.complete.obs",
      show_significance
    )
    
    meta_mat <- meta_res$correlation[ordered_features, , drop = FALSE]
    meta_signif <- meta_res$signif_matrix[ordered_features, , drop = FALSE]
    
    ht_meta <- create_cor_ht(
      meta_mat,
      meta_signif,
      "Environmental Variables",
      cluster_metadata_rows,
      cluster_metadata_columns
    )
    
    ht_list <- ht_meta
    
    if (verbose) message("✔ Metadata complete")
  }
  
  # ----------------------------
  # Reference heatmap
  # ----------------------------
  if (verbose) message("Generating reference heatmap...")
  
  ref_mat <- t(reference_layer[, ordered_features, drop = FALSE])
  
  ht_ref <- ComplexHeatmap::Heatmap(
    ref_mat,
    name = "Expression",
    col = expr_col_fun,
    cluster_rows = cluster_reference_columns,
    cluster_columns = cluster_reference_rows,
    show_row_names = show_column_names,
    show_column_names = show_row_names,
    row_names_gp = grid::gpar(fontsize = font_size),
    column_names_gp = grid::gpar(fontsize = font_size),
    column_title = "Transcriptomics",
    width = grid::unit(min(nrow(reference_layer) * 3, 250), "mm")
  )
  
  ht_list <- if (is.null(ht_list)) ht_ref else ht_list + ht_ref
  
  if (verbose) message("✔ Reference heatmap complete")
  
  # ----------------------------
  # Comparison layers
  # ----------------------------
  for (nm in names(comparison_layers)) {
    
    if (verbose)
      message("Computing correlations with ", nm, "...")
    
    res <- calculate_correlations(
      reference_layer,
      comparison_layers[[nm]],
      method,
      adjust_method,
      use,
      show_significance
    )
    
    cor_mat <- res$correlation[ordered_features, , drop = FALSE]
    signif_mat <- res$signif_matrix[ordered_features, , drop = FALSE]
    
    ht_layer <- create_cor_ht(
      cor_mat,
      signif_mat,
      nm,
      cluster_comparison_rows,
      cluster_comparison_columns
    )
    
    ht_list <- ht_list + ht_layer
    
    if (verbose)
      message("✔ ", nm, " complete")
  }
  
  if (verbose)
    message("Drawing final integrated heatmap...")
  
  ComplexHeatmap::draw(
    ht_list,
    merge_legends = TRUE,
    heatmap_legend_side = "right"
  )
  
  if (verbose) {
    message("✔ OmniCorr integration complete")
    message("==================================================")
  }
  
  invisible(ht_list)
}