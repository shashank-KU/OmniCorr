# OmniCorr

<img width="6000" height="4200" alt="Omnicorr (1)" src="https://github.com/user-attachments/assets/17cd5615-122c-4248-92b0-ec4ccde4cf3f" />

<!-- badges: start -->

[![Publication](https://img.shields.io/badge/Publication-Bioinformatics%20Advances-blue)](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbag057/8488725)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioadv%2Fvbag057-green)](https://doi.org/10.1093/bioadv/vbag057)
[![R](https://img.shields.io/badge/R-%3E%3D%204.6.0-276DC3)](https://cran.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GitHub last commit](https://img.shields.io/github/last-commit/shashank-KU/OmniCorr)](https://github.com/shashank-KU/OmniCorr)
[![GitHub issues](https://img.shields.io/github/issues/shashank-KU/OmniCorr)](https://github.com/shashank-KU/OmniCorr/issues)

<!-- badges: end -->

OmniCorr is an R package for correlation-based integration and visualization of multi-omics datasets, including transcriptomics, metatranscriptomics, metagenomics, metaproteomics, and associated sample metadata. The package facilitates feature-level comparison across omics layers and generates aligned heatmap visualizations highlighting putative cross-omics associations.

OmniCorr focuses on:

1. Computing pairwise correlations between independently processed omics layers.

2. Organizing correlation results in a reproducible data structure.

3. Visualizing correlations using aligned, hierarchical heatmaps.

OmniCorr is **not a preprocessing or network inference pipeline**. Instead, it assumes that each omics layer has already been processed using appropriate domain-specific workflows, such as WGCNA for transcriptomics. This design allows OmniCorr to be used with different preprocessing and upstream analysis workflows.

With appropriate modification of the input datasets, OmniCorr can be adapted to different omics data types. Users can also customize the heatmap color schemes and other visualization settings.

## Scope and assumptions

Before using OmniCorr, users should ensure that:

* Each omics dataset has been independently preprocessed.

* Samples are matched across datasets.

* Input matrices are numeric with samples as rows and features as columns.

* Feature-level summaries, such as hub genes, pathways, or taxa, are provided as input when required by the analysis.

## Publication

OmniCorr has been peer-reviewed and published in *Bioinformatics Advances*.

**OmniCorr: An R package for visualizing putative host-microbiome interactions using multi-omics data.**
*Bioinformatics Advances* (2026).
https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbag057/8488725

If you use OmniCorr in your research, please cite this article.

## Installation

OmniCorr requires R (>= 4.6.0).

### Development version

The development version can be installed from GitHub:

```r
install.packages("pak")
pak::pak("shashank-KU/OmniCorr")
```

### Bioconductor

OmniCorr is being prepared for submission to Bioconductor. Following acceptance, the package will be available through Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("OmniCorr")
```

### Optional packages for examples

Some examples in this README use packages listed in `Suggests`, including `pheatmap`, `cowplot`, and `ggplot2`. These packages are not required for the core OmniCorr functions but can be installed if needed:

```r
install.packages(c("pheatmap", "cowplot", "ggplot2"))
```

### Load package and example data

```r
library(OmniCorr)

data(Metagenomics)
data(Transcriptomics)
data(Metatranscriptomics)
```

## Input data provenance

The example datasets included in OmniCorr are derived from typical upstream omics analysis workflows.

### Example: metatranscriptomics

Metatranscriptomics input can be generated from:

* `blockwiseModules()` from WGCNA for module detection.

* `chooseTopHubInEachModule()` for hub feature selection.

The resulting object contains module-level summaries and representative hub features.

Other omics layers, such as metagenomics, may originate from functional annotation or pathway profiling workflows. OmniCorr does not enforce a specific preprocessing strategy.

## Sample matching and ordering

All datasets **must contain identical samples in the same order**.

Check sample alignment:

```r
all(rownames(Transcriptomics) == rownames(Metagenomics))
```

If sample order differs, use the `check_sample_order()` helper function:

```r
df_list <- check_sample_order(Transcriptomics, Metagenomics)

Transcriptomics <- df_list[[1]]
Metagenomics <- df_list[[2]]
```

## Two ways to run OmniCorr

OmniCorr can be used in **two complementary ways**, depending on the level of control required:

| Workflow                     | Description                                                                                                                          | Best for                                                               |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------- |
| **Manual workflow**          | Users perform clustering, correlation calculation, and heatmap construction step-by-step using individual OmniCorr helper functions. | Maximum control, method development, or custom visualization pipelines |
| **Fully automated workflow** | A single function `run_omnicorr()` performs the complete integration workflow automatically.                                         | Fast analysis, reproducible workflows, and standard analyses           |

---

## Manual workflow

### Step 1: Feature clustering

The following steps use OmniCorr to integrate transcriptomics and metagenomics data and visualize the results. First, perform hierarchical clustering of the transcriptomics data:

```r
dendro <- hclust(
  d = as.dist(1 - WGCNA::bicor(Transcriptomics, maxPOutliers = 0.05)),
  method = "ward.D2"
)
```

### Step 2: Transcriptomics heatmap

```r
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

Use the `calculate_correlations()` function from OmniCorr to calculate Pearson correlations between the transcriptomics and metagenomics data:

```r
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

* a correlation matrix

* p-values

* adjusted p-values

* significance annotations

### Step 4: Correlation heatmaps

Create a color ramp for the heatmap using `colorRampPalette()` and `RColorBrewer`:

```r
heatmap_colors <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(6, "RdBu"))
)(51)
```

Generate heatmaps of the correlations using `pheatmap()` with the dendrogram from Step 1:

```r
hm_mg <- pheatmap::pheatmap(
  corr_mg$correlation,
  color = heatmap_colors,
  cluster_rows = dendro,
  cluster_cols = FALSE,
  display_numbers = corr_mg$signif_matrix,
  breaks = seq(-1, 1, length.out = 51),
  show_rownames = FALSE,
  legend = FALSE,
  main = "Metagenomics"
)

hm_mt <- pheatmap::pheatmap(
  corr_mt$correlation,
  color = heatmap_colors,
  cluster_rows = dendro,
  cluster_cols = FALSE,
  display_numbers = corr_mt$signif_matrix,
  breaks = seq(-1, 1, length.out = 51),
  show_rownames = TRUE,
  legend = TRUE,
  main = "Metatranscriptomics"
)
```

### Optional downstream visualization

The `cowplot` package can be used to combine the heatmaps into a single figure:

```r
cowplot::plot_grid(
  tx_heatmap$gtable,
  hm_mg$gtable,
  hm_mt$gtable,
  ncol = 3,
  align = "h",
  rel_widths = c(3, 1, 2)
)
```

![Omics Integration](https://user-images.githubusercontent.com/30895959/223124413-71981e48-a295-48cd-959a-8aec5e15d863.png)

### Step 5: Integrate an external heatmap

OmniCorr can correlate omics features with phenotypic or environmental variables.

```r
data(metadata)

all(row.names(metadata) == row.names(Transcriptomics))

corr_meta <- calculate_correlations(
  df1 = Transcriptomics,
  df2 = metadata,
  use = "pairwise.complete.obs",
  show_significance = "stars"
)
```

Generate the metadata correlation heatmap:

```r
hm_meta <- pheatmap::pheatmap(
  corr_meta$correlation,
  color = heatmap_colors,
  cluster_rows = dendro,
  cluster_cols = FALSE,
  display_numbers = corr_meta$signif_matrix,
  breaks = seq(-1, 1, length.out = 51),
  show_rownames = FALSE,
  legend = FALSE,
  main = "Environmental Variables"
)
```

The resulting heatmaps can be combined using `cowplot`:

```r
cowplot::plot_grid(
  hm_meta$gtable,
  tx_heatmap$gtable,
  hm_mg$gtable,
  hm_mt$gtable,
  ncol = 4,
  align = "h",
  rel_widths = c(1.5, 3.5, 1, 2)
)
```

![Rplot01](https://user-images.githubusercontent.com/30895959/223760509-8c3d8f8e-d232-4c0c-8832-9aa4c1ecf5d9.png)

## Fully automated workflow

OmniCorr provides a single high-level function, `run_omnicorr()`, that performs the complete integration workflow automatically.

The function internally performs:

* Sample alignment validation.

* Reference feature clustering, optional and enabled by default.

* Pairwise cross-omics correlation calculation.

* Multiple testing correction using false discovery rate (FDR) by default.

* Significance annotation using stars, p-values, or correlation values.

* Construction of aligned multi-panel heatmaps.

### Basic usage

Run the full integration with one command:

```r
run_omnicorr(
  reference_layer = Transcriptomics,
  reference_name = "Transcriptomics",
  comparison_layers = list(
    Metagenomics = Metagenomics,
    Metatranscriptomics = Metatranscriptomics
  ),
  metadata = metadata
)
```

<img width="1400" height="800" alt="Omnicorr" src="https://github.com/user-attachments/assets/0a5ca5b4-6110-4e6b-a7e0-784096a524d1" />

### Changing the reference omics layer

Any omics dataset can be used as the reference layer. For example, to use metatranscriptomics as the reference:

```r
run_omnicorr(
  reference_layer = Metatranscriptomics,
  reference_name = "Metatranscriptomics",
  comparison_layers = list(
    Transcriptomics = Transcriptomics,
    Metagenomics = Metagenomics
  ),
  metadata = metadata
)
```

### Running OmniCorr without metadata

Metadata are optional. OmniCorr can be run using only the reference and comparison omics layers:

```r
run_omnicorr(
  reference_layer = Transcriptomics,
  reference_name = "Transcriptomics",
  comparison_layers = list(
    Metagenomics = Metagenomics,
    Metatranscriptomics = Metatranscriptomics
  )
)
```

### Customizing heatmap titles

The names displayed in the heatmaps can be customized using `reference_name`, comparison layer names, and `metadata_name`:

```r
run_omnicorr(
  reference_layer = Transcriptomics,
  reference_name = "Host Gene Expression",
  comparison_layers = list(
    Microbiome = Metagenomics,
    Microbial_Activity = Metatranscriptomics
  ),
  metadata = metadata,
  metadata_name = "Environmental Variables"
)
```

### Changing the correlation method

OmniCorr supports Pearson, Spearman, and Kendall correlation.

#### Spearman correlation

```r
run_omnicorr(
  reference_layer = Transcriptomics,
  comparison_layers = list(
    Metagenomics = Metagenomics,
    Metatranscriptomics = Metatranscriptomics
  ),
  metadata = metadata,
  method = "spearman"
)
```

#### Kendall correlation

```r
run_omnicorr(
  reference_layer = Transcriptomics,
  comparison_layers = list(
    Metagenomics = Metagenomics,
    Metatranscriptomics = Metatranscriptomics
  ),
  metadata = metadata,
  method = "kendall"
)
```

### Adjusting label visibility

Users can control whether row or column labels are shown:

```r
run_omnicorr(
  reference_layer = Transcriptomics,
  comparison_layers = list(
    Metagenomics = Metagenomics,
    Metatranscriptomics = Metatranscriptomics
  ),
  metadata = metadata,
  show_row_names = FALSE,
  show_column_names = TRUE
)
```
