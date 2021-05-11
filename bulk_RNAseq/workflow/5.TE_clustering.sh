#!/bin/bash
# 5.TE_clustering.sh

### Transform TE expression to Z-score
mkdir ./data/Zscore_matrix
Rscript ./scripts/run_Zscore.R


### Perform k-means clustering of the Z-score differentially expressed TEs
mkdir ./data/Kmean_clusters
Rscript ./scripts/run_Kmean_clust.R

###############