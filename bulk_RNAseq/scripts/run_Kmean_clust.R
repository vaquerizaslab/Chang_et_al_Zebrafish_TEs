#!/usr/bin/env Rscript
#run_Kmean_clust.R
# Make K-mean clustering for 7 clusters. 
# Order K-mean clustering to match expression from early to late.  
# Save K-mean clustering.  
# Plot matrices for each clustering.  

suppressMessages(library(tictoc))

### Load Zscore matrix
tic("Load Z-score matrix")
Zscore_DIR <- "./data/Zscore_matrix"
Zscore_TEl_defrag_DE_selfExp <- readRDS(file.path(Zscore_DIR,"Zscore_TEl_defrag_DE_selfExp.min10.rds"))
toc()

### K-mean clustering
tic("K-mean clustering")
Kmean_clust <- kmeans(Zscore_TEl_defrag_DE_selfExp, centers=as.numeric(7), iter.max=500, nstart=50, algorithm="Lloyd")
toc()

### Save the K-mean object
tic("Save the K-mean object")
Kmean_clust_DIR <- "./data/Kmean_clusters"
saveRDS(Kmean_clust, file.path(Kmean_clust_DIR,"Kmean_DE_TEl_defrag_min10_SelfExp.rds"))
toc()

### Re-name clusters
# Cluster ordering functions
mean_expression_per_cluster <- function(mat_Zscore, cluster_list, verbose=FALSE){
	mat_Zscore_meanC <- data.frame()
	for (i in seq(1,max(cluster_list))) {
		if(verbose==TRUE){
			writeLines(paste0("Mean of cluster: ",i," (",length(cluster_list[which(cluster_list==i)]),")"))
		}
		if(i == 1){
			mat_Zscore_meanC <- colMeans(mat_Zscore[which(rownames(mat_Zscore) %in% names(cluster_list[which(cluster_list==i)])),])
		} else{
			mat_Zscore_meanC <- rbind(mat_Zscore_meanC, colMeans(mat_Zscore[which(rownames(mat_Zscore) %in% names(cluster_list[which(cluster_list==i)])),]))
		}
	}
	rownames(mat_Zscore_meanC) <- paste0("Cluster_",seq(1,max(cluster_list)))
	return(mat_Zscore_meanC)
}
mean_by_rep <- function(matrix, verbose=FALSE){
	sample_names <- unique(gsub("_[rR]ep_[1-9]+","",colnames(matrix)))
	for (sample in sample_names) {
		if(verbose==TRUE){
			writeLines(paste0("Calculating mean of: ",sample))
		}
		if(sample == sample_names[1]){
			matrix_mean <- rowMeans(matrix[,grep(sample,colnames(matrix))])
		} else{
			matrix_mean <- cbind(matrix_mean, rowMeans(matrix[,grep(sample,colnames(matrix))]))
		}
	}
	colnames(matrix_mean) <- sample_names
	df <- as.data.frame(matrix_mean)
	return(df)
}
order_cluster_by_expression <- function(mat_Zscore, Kmean_c, verbose=FALSE){
	# Given a matrix with the expression values and the list of clusters,
	# re-name the cluster number so that they match the proper order of the
	# heatmap plot.
	#
	K_list_cluster <- Kmean_c$cluster
	# Get the mean of cluster expression and replicates
	if(verbose==TRUE){
		mat_Zscore_mean_ClustRep <- mean_by_rep(mean_expression_per_cluster(mat_Zscore, K_list_cluster, verbose=TRUE), verbose=TRUE)
	} else {
		mat_Zscore_mean_ClustRep <- mean_by_rep(mean_expression_per_cluster(mat_Zscore, K_list_cluster, verbose=FALSE), verbose=FALSE)
	}
	# Order based on the mean of cluster expression and replicates
	all_cluster <- rownames(mat_Zscore_mean_ClustRep)
	cluster_order <- c()
	for (j in seq(1,ncol(mat_Zscore_mean_ClustRep))) {
		max_col <- max(mat_Zscore_mean_ClustRep[,j])
		highest_exp_cluster <- rownames(mat_Zscore_mean_ClustRep)[which(mat_Zscore_mean_ClustRep[,j]==max_col)]
		if( !(highest_exp_cluster %in% cluster_order) ){
			cluster_order <- c(cluster_order, highest_exp_cluster)
	  }
	}
	# Append those cluster that did not have any stage with the highest expression}
	remaining_clusters <- all_cluster[!(all_cluster %in% cluster_order)]
	cluster_order <- c(cluster_order, remaining_clusters)
	# Rename cluster 1-2-3 in order to A-B-C
	K_list_cluster_rename <- K_list_cluster
	for (clust in cluster_order) {
		c <- gsub("Cluster_","",clust)
		if(verbose==TRUE){
			writeLines(c)
	  }
		K_list_cluster_rename[which(K_list_cluster_rename==c)] <- LETTERS[which(cluster_order==clust)]
	}
	#
	return(K_list_cluster_rename)
}
#

# Re-name (order) clusters.
# From naming "Cluster_1" to "Cluster_A" to avoid confusing pre and post ordered lists.
tic("Ordering and re-naming clusters")
Kmean_clust_rename <- order_cluster_by_expression(Zscore_TEl_defrag_DE_selfExp, Kmean_clust, verbose=FALSE)
toc()

### Save ordered and re-name K_lists
tic("Save the K-mean object re-named")
saveRDS(Kmean_clust_rename, file.path(Kmean_clust_DIR,"Kmean_DE_TEl_defrag_min10_SelfExp_reorder.rds"))
toc()

### Plot Z-score matrix
tic(Plot Z-score matrix)
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(tictoc))
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
plot_Zscore_heatmap <- function(mat_Zscore, Kmean_c, title="Title", crop_color_map_MinMax){
  N_clusters <- as.numeric(N_clusters)
  crop_color_map_MinMax <- as.numeric(crop_color_map_MinMax)
  #
  for (CLUST in LETTERS[1:N_clusters]) {
    if (CLUST == "A") {
      mat_Zscore_sorted <- mat_Zscore[match(names(Kmean_c[which(Kmean_c==CLUST)]), rownames(mat_Zscore)),]
    } else {
      mat_Zscore_sorted <- rbind(mat_Zscore_sorted,
                                 mat_Zscore[match(names(Kmean_c[which(Kmean_c==CLUST)]), rownames(mat_Zscore)),] )
    }
  }
  #
  annot_row_clusters <- data.frame(Clusters=paste0("Cluster_",Kmean_c), 
                                   row.names=names(Kmean_c))
  cluster_colors <- list(Clusters=colorRampPalette(brewer.pal(8,"Dark2"))(N_clusters))
  names(cluster_colors$Clusters) <- paste0("Cluster_",LETTERS[1:N_clusters])
  # Run pheatmap
  writeLines(paste0("Minimum and maximum values from the color scale map will be cropped by ",crop_color_map_MinMax))
  crop_MIN_MAX_val <- crop_color_map_MinMax
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  MIN_MAX <- max(abs(min(mat_Zscore_sorted)), max(mat_Zscore_sorted))-crop_MIN_MAX_val
  myBreaks <- c(seq(-MIN_MAX, 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(MIN_MAX/paletteLength, MIN_MAX, length.out=floor(paletteLength/2)))
  clust_tbl <- table(Kmean_c); names(clust_tbl) <- paste0("Cluster_",names(clust_tbl))
  #cum_sum_clusters <- cumsum(clust_tbl[match(visual_cluster_order,names(clust_tbl))])
  cum_sum_clusters <- cumsum(as.vector(table(Kmean_c)))
  # Include N in parenthesis in the title
  title <- paste0(title," (",nrow(mat_Zscore_sorted),")")
  heatmap_Zscore_cluster <- pheatmap(mat_Zscore_sorted, color=myColor, breaks=myBreaks, cluster_cols=FALSE, cluster_rows=FALSE, show_rownames=FALSE,
           annotation_row = annot_row_clusters, annotation_colors = cluster_colors, annotation_legend = TRUE,
           gaps_row=cum_sum_clusters,#gaps_col=c(90,95,100)
           main=title )
  #
  return(heatmap_Zscore_cluster)
}
#
heatmap_Zscore <- plot_Zscore_heatmap(Zscore_TEl_defrag_DE_selfExp, Kmean_clust_rename, title="DE self expressed TE loci", crop_color_map_MinMax=2)
# Save heatmap
save_pheatmap_pdf(heatmap_Zscore, "heatmap_Zscore_TEl_defrag_DE_selfExp.pdf", width=10, height=7)
toc()

################