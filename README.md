
# OmniCorr

<!-- badges: start -->
<!-- badges: end -->

The code provided in this repository can be used as a framework for exploring correlations between different omics datasets, including transcriptomics, metatranscriptomics, metagenomics,metaproteomics, and others. By utilizing packages such as WGCNA, pheatmap, and cowplot, it provides a flexible and customizable way to generate integrated heatmaps and visualize correlations between datasets.

With some modifications to the input datasets, this code can be adapted to work with different omics data types. For example, in Step 3, users can substitute in their own data frames to calculate correlations between different datasets. Additionally, the color schemes used in the heatmaps can be adjusted to best suit the user's data.

Overall, this code provides a powerful tool for exploring correlations between different omics datasets, and can be adapted to work with a wide range of data types.

## Installation

Install R (version >= 3.6.0) and RStudio (optional)

### Dependencies

check to see if packages are installed. Install them if they are not, then load them into the R session or follow their official installation steps

```r
install_pak <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", "WGCNA", "pheatmap", "RColorBrewer", "cowplot", "devtools")
install_pak(packages)
```

You can install the development version of OmniCorr:

``` r
devtools::install_github("shashank-KU/OmniCorr")
```

### Load libraries
``` r
library(OmniCorr)
```

### Load example datasets
``` r
data(Metagenomics)
data(Transcriptomics)
data(Metatranscriptomics)
```

The input data for this package comes from the output of various analysis steps performed on different types of omics data such as transcriptomics, metatranscriptomics, metagenomics, proteomics, etc. For example, in the case of metatranscriptomics, the input data is generated from the output of the `blockwiseModules` function of the `WGCNA` package, which performs a network analysis of gene expression data to identify gene modules. The `chooseTopHubInEachModule` function is then used to identify the top hub genes in each module. The resulting data is stored in a list called Metatranscriptomics, which contains information about the identified modules and top hub genes.

Similarly, for other types of omics data, different analysis steps are performed to generate input data for the package. For example, in the case of metagenomics, the input data may be generated from the output of functional annotation and pathway analysis tools. The specific analysis steps and tools used may vary depending on the type of omics data being analyzed.


***Note:*** Before using the OmniCorr package, it is important to ensure that the row names of the transcriptomics and metagenomics data frames are matching. 
You can check this by running the command ```table(rownames(Transcriptomics) == rownames(Metagenomics))```. If the result is not all TRUE, then you will need to modify your data frames to ensure that the row names match.

If samples is not is same order, use the `CheckSampleOrder` function below to match and reorder the samples. This function takes two data frames, df1 and df2, as inputs, and checks if their rownames match and are in the same order. If the rownames do not match, it reorders the rows of both data frames to match, and then checks again if the rownames match. If they do not match after reordering, it throws an error message. If the rownames match, it returns a list with both data frames in the same order and with matching rownames.

```r
df_list <- CheckSampleOrder(Transcriptomics, Metagenomics)
Transcriptomics <- df_list[[1]]
Metagenomics <- df_list[[2]]
```
## Steps

The following steps use OmniCorr to integrate transcriptomics and metagenomics data and visualize the results:

### Step 1: Perform hierarchical clustering of transcriptomics data using `WGCNA`

Load the WGCNA package in R
Calculate the bicorrelation matrix of transcriptomics data with outlier removal using the `bicor()` function
Convert the matrix to a distance matrix using `as.dist()`
Perform hierarchical clustering on the distance matrix using `hclust()` with `ward.D2` method

``` r
dendro <- hclust(as.dist(1 - WGCNA::bicor(Transcriptomics, maxPOutliers = 0.05)), method = "ward.D2")
```
### Step 2: Generate a heatmap of transcriptomics data with dendrogram

Transpose the transcriptomics data using `t()`
Create a data frame from the transposed data using `data.frame()`
Generate a heatmap of the data using `pheatmap()` with the dendrogram from Step 1 and no column tree

``` r
result2 <- pheatmap::pheatmap(t(Transcriptomics), 
                   cluster_rows = dendro, 
                   cluster_cols = F, 
                   show_rownames = F, 
                   main = paste("Transcriptomics"))
```

### Step 3: Calculate correlations between omics data

Use the `calculate_correlations()` function from `OmniCorr` to calculate `Pearson` correlations between the transposed transcriptomics data and the metagenomics data

``` r
result3 <- calculate_correlations(df1 = Transcriptomics, 
                                 df2 = Metagenomics,
                                 show_significance = "stars" # Possible other values are "p_value" or "correlation"
                                 )
```
Use the `calculate_correlations()` function from `OmniCorr` to calculate `Pearson` correlations between the transposed transcriptomics data and the metatranscriptomics data
``` r
result3.1 <- calculate_correlations(df1 = Transcriptomics, 
                                 df2 = Metatranscriptomics,
                                 show_significance = "stars" # Possible other values are "p_value" or "correlation"
                                 )
```

### Step 4: Generate a heatmap of the correlations between transcriptomics with metagenomics and metatranscriptomics data

Create a color ramp for the heatmap using `colorRampPalette()`
Generate a heatmap of the correlations using `pheatmap()` with the dendrogram from [Step 1](https://github.com/shashank-KU/OmniCorr#step-1-perform-hierarchical-clustering-of-transcriptomics-data-using-wgcna), the color ramp, and no row tree
Add significant correlations to the heatmap using the display_numbers parameter from `pheatmap()`

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

Use the `plot_grid()` function from the cowplot package to combine the two heatmaps into a single figure with two columns
Set the relative widths of the two columns using the `rel_widths` parameter
Adjust the margins of the plot using `ggplot2::theme()` with the `plot.margin` parameter


``` r
cowplot::plot_grid(result2$gtable, 
                   result4$gtable,
                   result4.1$gtable,
                   ncol = 3,  align = 'h',
                   rel_widths = c(3, 1, 2)) + 
                   ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"))
```

![Omics Integration](https://user-images.githubusercontent.com/30895959/223124413-71981e48-a295-48cd-959a-8aec5e15d863.png)



### Step 6: Integrate the external heatmap

Follow step 3, Use the `calculate_correlations()` function from `OmniCorr` to calculate `Pearson` correlations between the transposed transcriptomics data and the external variables. The metadata used in this package is an example data frame that contains information about the samples used in the study. It includes 14 variables, such as body weight, length, sex, gut condition, bleeding, fat and liver scores, and hepatosomatic and cardiosomatic indices.

``` r
data(metadata)
all(row.names(metadata) == row.names(Transcriptomics))

result3.2 <- calculate_correlations(df1 = Transcriptomics, 
                                    df2 = metadata, 
                                    use = "pairwise.complete.obs", # default is "all.obs"
                                    show_significance = "stars" # Possible other values are "p_value" or "correlation"
                                    )
   
```
``` r
heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)

result4.2 <- pheatmap::pheatmap(result3.2$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,
                   cluster_rows = dendro,
                   #cutree_rows = row_cut,
                   display_numbers = result3.2$signif_matrix, 
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   show_rownames = F, legend = F,
                   labels_row = paste0(rownames(result3.2$correlation)),
                   labels_col = paste0(colnames(result3.2$correlation)),
                   main = paste("Environmental Variables"))
```
``` r
cowplot::plot_grid(result4.2$gtable, # External heatmap correlation
                   result2$gtable,  # Transcriptomics heatmap correlation
                   result4$gtable,  # Metagenomics heatmap correlation
                   result4.1$gtable,  # Metatranscriptomics heatmap correlation
                   ncol = 4, # change this based upon the number of omics heatmap
                   align = 'h',
                   rel_widths = c(1.5, 3.5, 1, 2)) + 
                   ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"))
```
![Rplot01](https://user-images.githubusercontent.com/30895959/223760509-8c3d8f8e-d232-4c0c-8832-9aa4c1ecf5d9.png)
