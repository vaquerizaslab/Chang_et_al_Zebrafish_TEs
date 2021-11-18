Code for the heatmap (Fig. 4A)
================

**\#this is the code for the analysis of single cell sequencing in the
manuscript entitled “A genomic portrait of zebrafish transposable
elements and their spatiotemporal embryonic expression”**<br> \#\#Code
for the heatmap<br> \#\#grep 34 TE markers identified from Seurat
(blastula/gastrula/segmentation\*\* batchcorrected\_AE.txt)

``` r
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
df <- read.delim("34_DETEs_3stages_AE.txt", sep = "\t", stringsAsFactors=F, header = TRUE, row.names=1)
df.numeric.2 <- df+1

#z-score
x<-apply(df.numeric.2,1,mean)
sd<-apply(df.numeric.2,1,sd)
df.zscore<-(df.numeric.2-x)/sd

df.zscore1 = df.zscore[, 1:13]
df.zscore2 = df.zscore[, 14:52]
df.zscore3 = df.zscore[, 53:88]

# K-means clustering
km = kmeans(df.zscore, centers = 2)$cluster
split <- paste0("Cluster\n", km)

# Heirarchical clustering
dist_mat <- dist(df.zscore, method='euclidean')
hclust_ward <- hclust(dist_mat, method = 'ward.D2')
cut_ward <- cutree(hclust_ward, k=2)
hsplit <- paste0("Cluster\n", cut_ward)
plot(hclust_ward)

Heatmap(df.zscore1,
        cluster_rows = TRUE,
        row_names_side ="left",
        split=hsplit,
        row_order = hclust_ward$order,
        row_names_gp = gpar(fontsize = 15),
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        col = colorRamp2(c(-1, 0, 1), c("skyblue3", "white", "lightgoldenrod1"))) +
  Heatmap(df.zscore2,
          cluster_rows = FALSE,
          row_order = hclust_ward$order,
          show_row_names = FALSE,
          cluster_columns = FALSE,
          show_column_dend = FALSE,
          col = colorRamp2(c(-1, 0, 1), c("skyblue3", "white", "lightgoldenrod1"))) +
  Heatmap(df.zscore3,
          cluster_rows = FALSE,
          row_order = hclust_ward$order,
          show_row_names = FALSE,
          cluster_columns = FALSE,
          show_column_dend = FALSE,
          col = colorRamp2(c(-1, 0, 1), c("skyblue3", "white", "lightgoldenrod1")))

```
