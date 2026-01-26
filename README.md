
# OmniCorr

<!-- badges: start -->
<!-- badges: end -->

OmniCorr is an R package–style framework for correlation-based integration and visualization of multi-omics datasets, including transcriptomics, metatranscriptomics, metagenomics, metaproteomics, and associated sample metadata. The framework is designed to facilitate interpretable, feature-level comparison across omics layers and to generate aligned heatmap visualizations highlighting putative cross-omics associations.

OmniCorr focuses on:

1. Computing pairwise correlations between independently processed omics layers

2. Organizing correlation results in a reproducible data structure

3. Visualizing correlations using aligned, hierarchical heatmaps

OmniCorr is **not a preprocessing or network inference pipeline**. Instead, it assumes that each omics layer has already been processed using appropriate domain-specific workflows (e.g. WGCNA for transcriptomics). This design ensures flexibility and interoperability with existing bioinformatics pipelines.
With some modifications to the input datasets, this code can be adapted to work with different omics data types. For example, in Step 3, users can substitute in their own data frames to calculate correlations between different datasets. Additionally, the color schemes used in the heatmaps can be adjusted to best suit the user's data.


## Scope and assumptions
Before using OmniCorr, users should ensure that:

- Each omics dataset has been independently preprocessed

- Samples are matched across datasets

- Input matrices are numeric with samples as rows and features as columns

- Feature-level summaries (e.g. hub genes, pathways, taxa) are provided as input

## Installation

Install R (version >= 3.6.0) and RStudio (optional)

### Dependencies

Some dependencies are available via Bioconductor:
```r
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("impute", "preprocessCore", "GO.db"))
```
Install required CRAN packages:
```r
install_pak <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}
packages <- c("ggplot2", "WGCNA", "pheatmap", "RColorBrewer", "cowplot", "devtools")
install_pak(packages)
```

### Install OmniCorr (development version)
``` r
devtools::install_github("shashank-KU/OmniCorr")
```

### Alternative installation (Linux, without devtools)

```r
wget https://github.com/shashank-KU/OmniCorr/archive/refs/heads/master.zip
unzip master.zip
mv OmniCorr-master/ OmniCorr
R CMD INSTALL OmniCorr
```

### Load package and example data
``` r
library(OmniCorr)
data(Metagenomics)
data(Transcriptomics)
data(Metatranscriptomics)
```
## Input data provenance
The example datasets included in OmniCorr are derived from typical upstream omics analysis workflows.

### Example: metatranscriptomics
Metatranscriptomics input is generated from:

- `blockwiseModules()` (WGCNA) for module detection

- `chooseTopHubInEachModule()` for hub feature selection
The resulting object contains module-level summaries and representative hub features.
Other omics layers (e.g. metagenomics) may originate from functional annotation or pathway profiling pipelines. OmniCorr does not enforce a specific preprocessing strategy.


## Sample matching and ordering (required)
All datasets **must contain identical samples in the same order**.

Check sample alignment:
```table(rownames(Transcriptomics) == rownames(Metagenomics))```

If sample order differs, use the helper function:
```r
df_list <- CheckSampleOrder(Transcriptomics, Metagenomics)
Transcriptomics <- df_list[[1]]
Metagenomics <- df_list[[2]]
```
## Core workflow
### Step 1: Feature clustering (reference layer)

The following steps use OmniCorr to integrate transcriptomics and metagenomics data and visualize the results. First perform hierarchical clustering of transcriptomics

``` r
dendro <- hclust(as.dist(1 - WGCNA::bicor(Transcriptomics, maxPOutliers = 0.05)), method = "ward.D2")
```

### Step 2: Transcriptomics heatmap

``` r
tx_heatmap <- pheatmap::pheatmap(
  t(Transcriptomics),
  cluster_rows = dendro,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Transcriptomics"
)
```

This heatmap provides the structural backbone for aligning downstream correlation heatmaps.

### Step 3: Cross-omics correlation analysis

Use the `calculate_correlations()` function from `OmniCorr` to calculate `Pearson` correlations between the transposed transcriptomics data and the metagenomics data

``` r
corr_mg <- calculate_correlations(
  df1 = Transcriptomics,
  df2 = Metagenomics,
  show_significance = "stars"
)
corr_mt <- calculate_correlations(
  df1 = Transcriptomics,
  df2 = Metatranscriptomics,
  show_significance = "stars"
)
```

Each call returns:

- correlation matrix

- p-values

- significance annotations

### Step 4: Correlation heatmaps

Create a color ramp for the heatmap using `colorRampPalette()`
Generate a heatmap of the correlations using `pheatmap()` with the dendrogram from [Step 1](https://github.com/shashank-KU/OmniCorr#step-1-perform-hierarchical-clustering-of-transcriptomics-data-using-wgcna), the color ramp, and no row tree
Add significant correlations to the heatmap using the display_numbers parameter from `pheatmap()`

``` r
heatmap_colors <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(6, "RdBu"))
)(51)

hm_mg <- pheatmap::pheatmap(
  corr_mg$correlation,
  color = heatmap_colors,
  cluster_rows = dendro,
  display_numbers = corr_mg$signif_matrix,
  breaks = seq(-1, 1, length.out = 51),
  show_rownames = FALSE,
  legend = FALSE,
  main = "Metagenomics")

hm_mt <- pheatmap::pheatmap(
  corr_mt$correlation,
  color = heatmap_colors,
  cluster_rows = dendro,
  display_numbers = corr_mt$signif_matrix,
  breaks = seq(-1, 1, length.out = 51),
  show_rownames = TRUE,
  legend = TRUE,
  main = "Metatranscriptomics")
```

### Optional downstream visualization (recommended, not required)

Use the `plot_grid()` function from the cowplot package to combine the two heatmaps into a single figure with two columns
Set the relative widths of the two columns using the `rel_widths` parameter
Adjust the margins of the plot using `ggplot2::theme()` with the `plot.margin` parameter


``` r
cowplot::plot_grid(
  tx_heatmap$gtable,
  hm_mg$gtable,
  hm_mt$gtable,
  ncol = 3,
  align = "h",
  rel_widths = c(3, 1, 2)
) +
ggplot2::theme(
  plot.margin = ggplot2::unit(c(1,1,1,1), "cm")
)

```

![Omics Integration](https://user-images.githubusercontent.com/30895959/223124413-71981e48-a295-48cd-959a-8aec5e15d863.png)



### Step 5: Integrate the external heatmap

OmniCorr can correlate omics features with phenotypic or environmental variables.

``` r
data(metadata)
all(row.names(metadata) == row.names(Transcriptomics))

corr_meta <- calculate_correlations(
  df1 = Transcriptomics,
  df2 = metadata,
  use = "pairwise.complete.obs", # default is "all.obs"
  show_significance = "stars" # Possible other values are "p_value" or "correlation"
)
```
``` r
hm_meta <- pheatmap::pheatmap(
  corr_meta$correlation,
  color = heatmap_colors,
  cluster_rows = dendro,
  display_numbers = corr_meta$signif_matrix,
  breaks = seq(-1, 1, length.out = 51),
  show_rownames = FALSE,
  legend = FALSE,
  main = "Environmental Variables")
```
``` r
cowplot::plot_grid(
  hm_meta$gtable, # External heatmap correlation
  tx_heatmap$gtable, # Transcriptomics heatmap correlation
  hm_mg$gtable, # Metagenomics heatmap correlation
  hm_mt$gtable, # Metatranscriptomics heatmap correlation
  ncol = 4,
  align = "h",
  rel_widths = c(1.5, 3.5, 1, 2))
```
![Rplot01](https://user-images.githubusercontent.com/30895959/223760509-8c3d8f8e-d232-4c0c-8832-9aa4c1ecf5d9.png)


## Optional downstream visualization: Circos (chord) diagram

In addition to aligned heatmaps, OmniCorr provides an optional Circos (chord) visualization to summarize statistically significant cross-omics associations in a compact network-style layout. This visualization is intended for exploratory interpretation and is not part of the core OmniCorr pipeline.

The Circos plot displays:

- Omics features (e.g. modules, hub genes, taxa) grouped by data type

- Optional sample-level metadata variables

- Significant cross-omics correlations passing a user-defined FDR threshold

- Link color and width representing correlation direction and magnitude

### Key characteristics

- Supports any number of omics layers

- Uses pairwise complete observations

- Applies FDR correction across all tested associations

- Retains metadata variables regardless of degree

- Designed for summary visualization, not statistical inference

### Example usage
``` r
omics <- list(
  Transcriptomics     = Transcriptomics,
  Metagenomics        = Metagenomics,
  Metatranscriptomics = Metatranscriptomics
)

plot_omnicorr_circos(
  omics_list = omics,
  metadata   = metadata,
  fdr_cutoff = 0.0001
)
```
