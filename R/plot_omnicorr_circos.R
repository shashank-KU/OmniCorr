#' Circos visualization of cross-omics correlations using OmniCorr
#'
#' @description
#' Generates an optional Circos (chord) diagram to visualize statistically
#' significant correlations between multiple omics layers and optional
#' numeric sample metadata.
#'
#' @param omics_list Named list of omics matrices (rows = samples, cols = features).
#' @param metadata Optional numeric sample-level metadata.
#' @param method Correlation method ("spearman", "pearson", "kendall").
#' @param fdr_cutoff Adjusted p-value cutoff.
#' @param max_abs_cor Maximum absolute correlation used for color scaling.
#' @param min_degree Minimum number of edges required for omics features.
#'
#' @return Invisibly returns a data frame of retained edges.
#'
#' @export
plot_omnicorr_circos <- function(
    omics_list,
    metadata = NULL,
    method = "spearman",
    fdr_cutoff = 0.05,
    max_abs_cor = 0.5,
    min_degree = 2
) {
  
  stopifnot(is.list(omics_list), length(omics_list) >= 2)
  
  ## ===========================
  ## Sample consistency checks
  ## ===========================
  if (any(vapply(omics_list, function(x) is.null(rownames(x)), logical(1)))) {
    stop("All omics matrices must have sample IDs as rownames.")
  }
  
  n_samples <- vapply(omics_list, nrow, integer(1))
  if (length(unique(n_samples)) != 1) {
    stop("All omics matrices must have the same number of samples.")
  }
  
  ref_samples <- rownames(omics_list[[1]])
  for (nm in names(omics_list)) {
    if (!identical(rownames(omics_list[[nm]]), ref_samples)) {
      stop(
        sprintf(
          "Sample mismatch detected in omics layer '%s'.\nUse CheckSampleOrder() before calling plot_omnicorr_circos().",
          nm
        )
      )
    }
  }
  
  if (!is.null(metadata)) {
    if (is.null(rownames(metadata))) {
      stop("Metadata must have sample IDs as rownames.")
    }
    if (!identical(rownames(metadata), ref_samples)) {
      stop(
        "Sample mismatch between metadata and omics data.\nUse CheckSampleOrder() before calling plot_omnicorr_circos()."
      )
    }
  }
  
  ## ===========================
  ## Safe correlation helpers
  ## ===========================
  safe_cor <- function(x, y, method) {
    ok <- stats::complete.cases(x, y)
    if (sum(ok) < 3) return(NA_real_)
    stats::cor(x[ok], y[ok], method = method)
  }
  
  safe_cor_test <- function(x, y, method) {
    ok <- stats::complete.cases(x, y)
    if (sum(ok) < 3) return(NA_real_)
    suppressWarnings(stats::cor.test(x[ok], y[ok], method = method)$p.value)
  }
  
  ## ===========================
  ## Helper: compute edges
  ## ===========================
  compute_edges <- function(df1, df2, name1, name2) {
    
    expand.grid(
      f1 = colnames(df1),
      f2 = colnames(df2),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        rho  = safe_cor(df1[[f1]], df2[[f2]], method),
        pval = safe_cor_test(df1[[f1]], df2[[f2]], method)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        fdr  = stats::p.adjust(pval, method = "fdr"),
        from = paste(name1, f1, sep = "::"),
        to   = paste(name2, f2, sep = "::")
      ) %>%
      dplyr::filter(!is.na(rho), fdr <= fdr_cutoff) %>%
      dplyr::select(from, to, rho)
  }
  
  ## ===========================
  ## Omics–omics correlations
  ## ===========================
  edges <- list()
  omic_names <- names(omics_list)
  
  for (i in seq_len(length(omics_list) - 1)) {
    for (j in (i + 1):length(omics_list)) {
      edges[[paste(i, j, sep = "_")]] <-
        compute_edges(
          omics_list[[i]],
          omics_list[[j]],
          omic_names[i],
          omic_names[j]
        )
    }
  }
  
  edges <- dplyr::bind_rows(edges)
  
  ## ===========================
  ## Metadata correlations
  ## ===========================
  if (!is.null(metadata)) {
    
    metadata <- metadata[, vapply(metadata, is.numeric, logical(1)), drop = FALSE]
    
    edges_meta <- lapply(names(omics_list), function(nm) {
      
      expand.grid(
        f1 = colnames(omics_list[[nm]]),
        f2 = colnames(metadata),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          rho  = safe_cor(omics_list[[nm]][[f1]], metadata[[f2]], method),
          pval = safe_cor_test(omics_list[[nm]][[f1]], metadata[[f2]], method)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          fdr  = stats::p.adjust(pval, method = "fdr"),
          from = paste(nm, f1, sep = "::"),
          to   = paste("Metadata", f2, sep = "::")
        ) %>%
        dplyr::filter(!is.na(rho), fdr <= fdr_cutoff) %>%
        dplyr::select(from, to, rho)
    }) %>% dplyr::bind_rows()
    
    edges <- dplyr::bind_rows(edges, edges_meta)
  }
  
  ## ===========================
  ## Truncate correlations
  ## ===========================
  edges <- edges %>%
    dplyr::mutate(
      rho_trunc = pmax(pmin(rho, max_abs_cor), -max_abs_cor)
    )
  
  ## ===========================
  ## Degree filtering
  ## ===========================
  node_degree <- edges %>%
    dplyr::select(from, to) %>%
    tidyr::pivot_longer(cols = c(from, to), values_to = "node") %>%
    dplyr::count(node)
  
  metadata_nodes <- node_degree$node[grepl("^Metadata::", node_degree$node)]
  omic_nodes <- node_degree %>%
    dplyr::filter(!grepl("^Metadata::", node) & n >= min_degree) %>%
    dplyr::pull(node)
  
  keep_nodes <- union(metadata_nodes, omic_nodes)
  
  edges <- edges %>%
    dplyr::filter(from %in% keep_nodes & to %in% keep_nodes)
  
  ## ===========================
  ## Circos plot
  ## ===========================
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    gap.degree = 1,
    track.margin = c(0.005, 0.005),
    canvas.xlim = c(-1.15, 1.15),
    canvas.ylim = c(-1.15, 1.15),
    cell.padding = c(0.001, 0, 0.001, 0)
  )
  
  sectors <- unique(c(edges$from, edges$to))
  
  circlize::circos.initialize(
    factors = sectors,
    xlim = cbind(rep(0, length(sectors)), rep(1, length(sectors)))
  )
  
  group <- ifelse(grepl("^Metadata::", sectors), "Metadata", sub("::.*", "", sectors))
  group_cols <- stats::setNames(
    RColorBrewer::brewer.pal(max(3, length(unique(group))), "Set2")[seq_along(unique(group))],
    unique(group)
  )
  sector_cols <- stats::setNames(group_cols[group], sectors)
  
  col_fun <- circlize::colorRamp2(
    c(-max_abs_cor, 0, max_abs_cor),
    c("#4575B4", "#F7F7F7", "#D73027")
  )
  
  circlize::chordDiagram(
    edges[, c("from", "to")],
    grid.col = sector_cols,
    col = col_fun(edges$rho_trunc),
    transparency = 0.3,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.05),
    link.lwd = pmax(abs(edges$rho_trunc) * 4, 0.5),
    directional = 0
  )
  
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      label <- sub("^.*::", "", circlize::CELL_META$sector.index)
      circlize::circos.text(
        mean(circlize::CELL_META$xlim),
        circlize::CELL_META$ylim[1] + circlize::mm_y(2),
        label,
        facing = "clockwise",
        niceFacing = TRUE,
        cex = 0.6
      )
    },
    bg.border = NA
  )
  
  lgd_group <- ComplexHeatmap::Legend(
    labels = names(group_cols),
    legend_gp = grid::gpar(fill = group_cols),
    title = "Data type"
  )
  
  lgd_cor <- ComplexHeatmap::Legend(
    col_fun = col_fun,
    title = paste(method, "correlation")
  )
  
  ComplexHeatmap::draw(
    ComplexHeatmap::packLegend(lgd_group, lgd_cor),
    x = grid::unit(1, "npc") - grid::unit(1.5, "cm"),
    y = grid::unit(0.5, "npc"),
    just = c("right", "center")
  )
  
  invisible(edges)
}
