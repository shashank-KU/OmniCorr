#' Circos visualization of cross-omics correlations using OmniCorr
#'
#' @description
#' Generates an optional Circos (chord) diagram to visualize statistically
#' significant correlations between multiple omics layers (e.g. transcriptomics,
#' metagenomics, metatranscriptomics) and optional numeric sample metadata.
#'
#' This function is intended as a downstream visualization utility and is not
#' part of the core OmniCorr correlation or clustering pipeline.
#'
#' @param omics_list A named list of data frames or matrices, where each element
#'   represents one omics layer. Rows must correspond to samples and columns to
#'   features. All omics layers must share identical sample order.
#'
#' @param metadata (optional) A data frame or matrix of numeric sample-level
#'   metadata (rows = samples, columns = variables).
#'
#' @param method Correlation method. Default is \code{"spearman"}.
#'
#' @param fdr_cutoff Adjusted p-value threshold for retaining correlations.
#'   Default is \code{0.05}.
#'
#' @param max_abs_cor Maximum absolute correlation value used for color scaling.
#'   Default is \code{0.5}.
#'
#' @param min_degree Minimum number of significant connections required for
#'   omics features to be displayed. Metadata variables are always retained.
#'
#' @return Invisibly returns a data frame of retained edges used in the Circos plot.
#'
#' @examples
#' \dontrun{
#' omics <- list(
#'   Transcriptomics = Transcriptomics,
#'   Metagenomics = Metagenomics,
#'   Metatranscriptomics = Metatranscriptomics
#' )
#'
#' plot_omnicorr_circos(
#'   omics_list = omics,
#'   metadata = metadata
#' )
#' }
#'
#' @importFrom stats cor cor.test p.adjust complete.cases
#' @importFrom dplyr %>% mutate filter select count bind_rows rowwise ungroup pull
#' @importFrom tidyr pivot_longer
#' @importFrom circlize circos.clear circos.par circos.initialize circos.track chordDiagram colorRamp2 mm_y CELL_META
#' @importFrom ComplexHeatmap Legend packLegend draw
#' @importFrom grid unit gpar
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
    stats::cor.test(x[ok], y[ok], method = method)$p.value
  }
  
  ## ===========================
  ## Compute omics–omics edges
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
  ## Pairwise omics correlations
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
    
    metadata <- as.data.frame(metadata)
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
  ## Circos setup
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
  
  ## ===========================
  ## Groups & colors
  ## ===========================
  group <- ifelse(
    grepl("^Metadata::", sectors),
    "Metadata",
    sub("::.*", "", sectors)
  )
  
  group_levels <- unique(group)
  
  group_cols <- stats::setNames(
    RColorBrewer::brewer.pal(max(3, length(group_levels)), "Set2")[seq_along(group_levels)],
    group_levels
  )
  
  sector_cols <- stats::setNames(group_cols[group], sectors)
  
  col_fun <- circlize::colorRamp2(
    c(-max_abs_cor, 0, max_abs_cor),
    c("#4575B4", "#F7F7F7", "#D73027")
  )
  
  ## ===========================
  ## Draw Circos
  ## ===========================
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
  
  ## ===========================
  ## Feature labels
  ## ===========================
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      label <- sub("^.*::", "", circlize::CELL_META$sector.index)
      circlize::circos.text(
        x = mean(circlize::CELL_META$xlim),
        y = circlize::CELL_META$ylim[1] + circlize::mm_y(2),
        labels = label,
        facing = "clockwise",
        niceFacing = TRUE,
        cex = 0.6
      )
    },
    bg.border = NA
  )
  
  ## ===========================
  ## Legends (right side)
  ## ===========================
  lgd_group <- ComplexHeatmap::Legend(
    labels = names(group_cols),
    legend_gp = grid::gpar(fill = group_cols),
    title = "Data type"
  )
  
  lgd_cor <- ComplexHeatmap::Legend(
    col_fun = col_fun,
    title = paste(method, "correlation"),
    direction = "vertical"
  )
  
  ComplexHeatmap::draw(
    ComplexHeatmap::packLegend(lgd_group, lgd_cor),
    x = grid::unit(1, "npc") - grid::unit(1.5, "cm"),
    y = grid::unit(0.5, "npc"),
    just = c("right", "center")
  )
  
  invisible(edges)
}