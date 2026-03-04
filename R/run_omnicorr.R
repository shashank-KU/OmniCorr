#' Run OmniCorr Multi-Omics Integration
#'
#' @description
#' Integrates multiple omics layers by computing correlations between a
#' reference omics dataset and multiple comparison layers and visualizing the
#' results as aligned multi-panel heatmaps using \pkg{ComplexHeatmap}.
#'
#' The reference layer is displayed as an expression heatmap with optional
#' clustering of features and samples. Correlation heatmaps for comparison
#' layers and metadata are aligned to the same reference feature order to
#' ensure consistent interpretation across panels.
#'
#' @param reference_layer Numeric matrix/data.frame (samples x features).
#' @param comparison_layers Named list of numeric matrices/data.frames.
#' @param metadata Optional numeric matrix/data.frame (samples x variables).
#' #' @param reference_name Title for the reference heatmap.
#' @param comparison_names Optional titles for comparison heatmaps.
#' @param metadata_name Title for metadata heatmap.
#' @param method Correlation method: `"pearson"`, `"spearman"`, `"kendall"`.
#' @param adjust_method P-value correction method (`stats::p.adjust.methods`).
#' @param use Missing value handling strategy (`"all.obs"` or `"pairwise.complete.obs"`).
#' @param show_significance `"stars"`, `"p_value"`, or `"correlation"`.
#' @param cluster_reference_rows Cluster reference features.
#' @param cluster_reference_columns Cluster reference samples.
#' @param show_row_names Logical; display row labels.
#' @param show_column_names Logical; display column labels.
#' @param font_size Font size for row/column labels.
#' @param star_font_size Font size for cell annotations.
#' @param verbose Print progress messages.
#'
#' @return Invisibly returns a `ComplexHeatmap::HeatmapList`.
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar grid.text
#' @importFrom stats cor hclust as.dist
#' @importFrom grDevices colorRampPalette
#'
#' @export
run_omnicorr <- function(
    reference_layer,
    comparison_layers,
    metadata = NULL,
    reference_name = "Reference Layer",
    comparison_names = NULL,
    metadata_name = "Metadata",
    method = "pearson",
    adjust_method = "fdr",
    use = "all.obs",
    show_significance = "stars",
    cluster_reference_rows = TRUE,
    cluster_reference_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    font_size = 8,
    star_font_size = 8,
    verbose = TRUE
){
  
  required_pkgs <- c("ComplexHeatmap","circlize","RColorBrewer","grid")
  
  for(pkg in required_pkgs){
    if(!requireNamespace(pkg, quietly = TRUE)){
      stop(sprintf("Package '%s' must be installed.", pkg))
    }
  }
  
  method <- match.arg(method, c("pearson","spearman","kendall"))
  use <- match.arg(use, c("all.obs","pairwise.complete.obs"))
  show_significance <- match.arg(show_significance,
                                 c("stars","p_value","correlation"))
  adjust_method <- match.arg(adjust_method,
                             c("none", stats::p.adjust.methods))
  
  reference_layer <- as.matrix(reference_layer)
  
  if(is.null(rownames(reference_layer)))
    stop("reference_layer must contain rownames (sample IDs).")
  
  if(is.null(colnames(reference_layer)))
    stop("reference_layer must contain feature names.")
  
  if(!is.list(comparison_layers))
    stop("comparison_layers must be a named list.")
  
  if(is.null(names(comparison_layers)))
    stop("comparison_layers must have names.")
  
  if(is.null(comparison_names))
    comparison_names <- names(comparison_layers)
  
  if(length(comparison_names) != length(comparison_layers))
    stop("comparison_names length must match comparison_layers.")
  
  for(nm in names(comparison_layers)){
    layer <- as.matrix(comparison_layers[[nm]])
    
    if(!identical(rownames(layer), rownames(reference_layer)))
      stop(sprintf("Sample order mismatch in layer '%s'.", nm))
    
    comparison_layers[[nm]] <- layer
  }
  
  if(!is.null(metadata)){
    metadata <- as.matrix(metadata)
    
    if(!identical(rownames(metadata), rownames(reference_layer)))
      stop("Sample order mismatch between metadata and reference_layer.")
  }
  
  if(use == "all.obs"){
    if(anyNA(reference_layer) ||
       any(sapply(comparison_layers, anyNA)) ||
       (!is.null(metadata) && anyNA(metadata))){
      warning("Missing values detected. Switching to 'pairwise.complete.obs'.")
      use <- "pairwise.complete.obs"
    }
  }
  
  if(verbose){
    message("==================================================")
    message("OmniCorr Multi-Omics Integration")
    message("==================================================")
    message("Samples: ", nrow(reference_layer))
    message("Reference Features: ", ncol(reference_layer))
    message("Sample alignment verified")
  }
  
  feature_dendro <- FALSE
  ordered_features <- colnames(reference_layer)
  
  if(cluster_reference_rows){
    
    if(verbose) message("Clustering reference features...")
    
    cor_ref <- stats::cor(reference_layer, method = method, use = use)
    cor_ref[is.na(cor_ref)] <- 0
    
    feature_dendro <- stats::hclust(stats::as.dist(1 - cor_ref))
    ordered_features <- colnames(reference_layer)[feature_dendro$order]
  }
  
  sample_dendro <- FALSE
  
  if(cluster_reference_columns){
    
    if(verbose) message("Clustering samples...")
    
    cor_samples <- stats::cor(t(reference_layer), method = method, use = use)
    cor_samples[is.na(cor_samples)] <- 0
    
    sample_dendro <- stats::hclust(stats::as.dist(1 - cor_samples))
  }
  
  rdBu_colors <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(6,"RdBu"))
  )(51)
  
  cor_col_fun <- circlize::colorRamp2(seq(-1,1,length.out=51), rdBu_colors)
  
  expr_range <- range(reference_layer, na.rm = TRUE)
  
  expr_col_fun <- circlize::colorRamp2(
    seq(expr_range[1],expr_range[2],length.out=51),
    rdBu_colors
  )
  
  create_cor_ht <- function(mat, signif_mat, title_text){
    
    ComplexHeatmap::Heatmap(
      mat,
      name = paste0("Correlation_", title_text),
      col = cor_col_fun,
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = grid::gpar(fontsize = font_size),
      column_names_gp = grid::gpar(fontsize = font_size),
      column_title = title_text,
      
      cell_fun = function(j,i,x,y,width,height,fill){
        
        txt <- if(show_significance == "stars"){
          signif_mat[i,j]
        } else if(show_significance == "p_value"){
          sprintf("%.2g", signif_mat[i,j])
        } else {
          sprintf("%.2f", signif_mat[i,j])
        }
        
        grid::grid.text(txt,x,y,
                        gp = grid::gpar(fontsize = star_font_size))
      }
    )
  }
  
  ht_list <- list()
  
  if(!is.null(metadata)){
    
    if(verbose) message("Computing metadata correlations...")
    
    meta_res <- calculate_correlations(
      reference_layer,
      metadata,
      method,
      adjust_method,
      use,
      show_significance
    )
    
    meta_mat <- meta_res$correlation[ordered_features,,drop=FALSE]
    meta_sig <- meta_res$signif_matrix[ordered_features,,drop=FALSE]
    
    ht_list <- c(ht_list,
                 create_cor_ht(meta_mat, meta_sig, metadata_name))
  }
  
  if(verbose) message("Generating reference heatmap...")
  
  ref_mat <- t(reference_layer[,ordered_features,drop=FALSE])
  
  ht_ref <- ComplexHeatmap::Heatmap(
    ref_mat,
    name = paste0("Correlation_", reference_name),
    col = expr_col_fun,
    cluster_rows = feature_dendro,
    cluster_columns = sample_dendro,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = grid::gpar(fontsize = font_size),
    column_names_gp = grid::gpar(fontsize = font_size),
    column_title = reference_name
  )
  
  ht_list <- c(ht_list, ht_ref)
  
  for(i in seq_along(comparison_layers)){
    
    nm <- names(comparison_layers)[i]
    
    if(verbose)
      message("Computing correlations with ", nm, "...")
    
    res <- calculate_correlations(
      reference_layer,
      comparison_layers[[i]],
      method,
      adjust_method,
      use,
      show_significance
    )
    
    cor_mat <- res$correlation[ordered_features,,drop=FALSE]
    sig_mat <- res$signif_matrix[ordered_features,,drop=FALSE]
    
    ht_list <- c(ht_list,
                 create_cor_ht(cor_mat, sig_mat, comparison_names[i]))
  }
  
  if(verbose) message("Drawing final integrated heatmap...")
  
  ht <- Reduce(`+`, ht_list)
  
  ComplexHeatmap::draw(
    ht,
    merge_legends = TRUE,
    heatmap_legend_side = "right"
  )
  
  if(verbose){
    message("OmniCorr integration complete")
    message("==================================================")
  }
  
  invisible(ht)
  
}