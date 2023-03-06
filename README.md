
# OmicsIntegrator

<!-- badges: start -->
<!-- badges: end -->

The code provided in this repository can be used as a framework for exploring correlations between different omics datasets, including transcriptomics, metatranscriptomics, metagenomics,metaproteomics, and others. By utilizing packages such as WGCNA, pheatmap, and cowplot, it provides a flexible and customizable way to generate integrated heatmaps and visualize correlations between datasets.

With some modifications to the input datasets, this code can be adapted to work with different omics data types. For example, in Step 3, users can substitute in their own data frames to calculate correlations between different datasets. Additionally, the color schemes used in the heatmaps can be adjusted to best suit the user's data.

Overall, this code provides a powerful tool for exploring correlations between different omics datasets, and can be adapted to work with a wide range of data types.

## Installation

Install R (version >= 3.6.0) and RStudio (optional)
Install the devtools package in R using install.packages("devtools")

You can install the development version of OmicsIntegrator like so:

``` r
devtools::install_github("shashank-KU/OmicsIntegrator")
```
Load libraries
``` r
library(OmicsIntegrator)
library(pheatmap)
```

Load example datasets
``` r
data(Metagenomics)
data(Transcriptomics)
data(Metatranscriptomics)
```

## Steps

The following steps use OmicsIntegrator to integrate transcriptomics and metagenomics data and visualize the results:

### Step 1: Perform hierarchical clustering of transcriptomics data using WGCNA

Load the WGCNA package in R
Calculate the bicorrelation matrix of transcriptomics data with outlier removal using the bicor() function
Convert the matrix to a distance matrix using as.dist()
Perform hierarchical clustering on the distance matrix using hclust() with "ward.D2" method

``` r
dendro <- hclust(as.dist(1 - WGCNA::bicor(Transcriptomics, maxPOutliers = 0.05)), method = "ward.D2")
```

### Step 2: Generate a heatmap of transcriptomics data with dendrogram

Transpose the transcriptomics data using t()
Create a data frame from the transposed data using data.frame()
Generate a heatmap of the data using pheatmap() with the dendrogram from Step 1 and no column tree

``` r
Transcriptomics<-data.frame(t(Transcriptomics))
result2 <-pheatmap(Transcriptomics, 
                   cluster_rows = dendro, 
                   cluster_cols = F, 
                   show_rownames = F, 
                   main = paste("Transcriptomics"))
```

### Step 3: Calculate correlations between transcriptomics and metagenomics data using OmicsIntegrator

Use the calculate_correlations() function from OmicsIntegrator to calculate Pearson correlations between the transposed transcriptomics data and the metagenomics data

``` r
result3 <- calculate_correlations(df1 = t(Transcriptomics), 
                                 df2 = Metagenomics)
```
Use the calculate_correlations() function from OmicsIntegrator to calculate Pearson correlations between the transposed transcriptomics data and the metatranscriptomics data
``` r
result3.1 <- calculate_correlations(df1 = t(Transcriptomics), 
                                 df2 = Metatranscriptomics)
```

### Step 4: Generate a heatmap of the correlations between transcriptomics and metagenomics data

Create a color ramp for the heatmap using colorRampPalette()
Generate a heatmap of the correlations using pheatmap() with the dendrogram from Step 1, the color ramp, and no row tree
Add significant correlations to the heatmap using the display_numbers parameter from pheatmap()

``` r
heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)

result4 <- pheatmap::pheatmap(result3$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,
                   cluster_rows = dendro,
                   #cutree_rows = row_cut,
                   display_numbers = result3$signif_matrix, 
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   show_rownames = F, legend = F,
                   labels_row = paste0(rownames(result3$correlation)),
                   labels_col = paste0(colnames(result3$correlation)),
                   main = paste("Metagenomics"))

result4.1 <- pheatmap::pheatmap(result3.1$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,
                   cluster_rows = dendro,
                   #cutree_rows = row_cut,
                   display_numbers = result3.1$signif_matrix, 
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   show_rownames = T, legend = T,
                   labels_row = paste0(rownames(result3.1$correlation)),
                   labels_col = paste0(colnames(result3.1$correlation)),
                   main = paste("Metatranscriptomics"))
```

### Step 5: Combine the heatmap of transcriptomics data and the heatmap of correlations

Use the plot_grid() function from the cowplot package to combine the two heatmaps into a single figure with two columns
Set the relative widths of the two columns using the rel_widths parameter
Adjust the margins of the plot using ggplot2::theme() with the plot.margin parameter


``` r
cowplot::plot_grid(result2$gtable, 
                   result4$gtable,
                   result4.1$gtable,
                   ncol = 3,  align = 'h',
                   rel_widths = c(3, 1, 2)) + 
  ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"))
```

![Omics Integration](https://user-images.githubusercontent.com/30895959/223124413-71981e48-a295-48cd-959a-8aec5e15d863.png)





