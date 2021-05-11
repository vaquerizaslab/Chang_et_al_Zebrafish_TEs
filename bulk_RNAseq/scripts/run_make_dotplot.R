#!/usr/bin/env Rscript
#run_make_dotplot.R
# Perform TE enrichment analysis

suppressMessages(library(tictoc))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggplot2))

### Load TE annotations with extended 3' UTR data
tic("Loading TE annotation")
TE_GTF_DFRAG_path <- "./data/danRer11.TEtrans_uID.dfrag.classified.gtf"
TE_GTF_defrag <- import(TE_GTF_DFRAG_path)
toc(); writeLines("")

### Load ordered and re-name K-mean clustering
tic("Load ordered and re-name K-mean clustering")
Kmean_clust_DIR <- "./data/Kmean_clusters"
Kmean_clust_rename <- readRDS(file.path(Kmean_clust_DIR,"Kmean_DE_TEl_defrag_min10_SelfExp_reorder.rds"))
# Remove TEs that are not in TE_GTF_defrag from Kmean_clust_rename files
# (Does are TEs in unplaced scaffolds)
Kmean_clust_rename <- Kmean_clust_rename[which( names(Kmean_clust_rename) %in% TE_GTF_defrag$defrag_transcript_id )]
toc()

### Functions for enrichment analysis
get_contTable_vals <- function(TEtype, TE, CLUSTER, KMER_DF){
	Kmer_DF_tmp <- KMER_DF
	Kmer_DF_tmp$Clust_tmp <- as.character(Kmer_DF_tmp$Cluster)
	if(TEtype=="TEclass"){
		Kmer_DF_tmp$TE_tmp <- as.character(Kmer_DF_tmp$TEclass)
	} else if(TEtype=="TEfam"){
		Kmer_DF_tmp$TE_tmp <- as.character(Kmer_DF_tmp$TEfam)
	} else if(TEtype=="TEsubfam"){
		Kmer_DF_tmp$TE_tmp <- as.character(Kmer_DF_tmp$TEsubfam)
	} else{
		stop("TEtype argument must be \"TEclass\", \"TEfam\" or \"TEsubfam\".")
	}
	Kmer_DF_tmp$Clust_tmp[which(Kmer_DF_tmp$Clust_tmp!=CLUSTER)] <- paste0("Not_",CLUSTER)
	Kmer_DF_tmp$TE_tmp[which(Kmer_DF_tmp$TE_tmp!=TE)] <- paste0("Not_",TE)
	contTable <- as.data.frame(table(Kmer_DF_tmp$TE_tmp,Kmer_DF_tmp$Clust_tmp))
	#print(contTable)
	contTable_vals <- c(contTable$Freq[which( contTable$Var1==paste0(TE)        & contTable$Var2==paste0(CLUSTER) )],
						contTable$Freq[which( contTable$Var1==paste0(TE)        & contTable$Var2==paste0("Not_",CLUSTER) )],
						contTable$Freq[which( contTable$Var1==paste0("Not_",TE) & contTable$Var2==paste0(CLUSTER) )],
						contTable$Freq[which( contTable$Var1==paste0("Not_",TE) & contTable$Var2==paste0("Not_",CLUSTER) )] )
	names(contTable_vals) <- c("TE_C","TE_noC","noTE_C","noTE_noC")
	#
	return(contTable_vals)
}
fisher_test_cluster_TEloci <- function(Kmean_c, N_clusters, TE_type){
	# Return a data frame with a line for each TE class/cluster,
	# and the results of the fisher exact test (p-value and Odds Ratio),
	# and the p-adjusted value.
	Kmer_vector <- Kmean_c
	Kmer_DF <- data.frame(Cluster=Kmer_vector)
	Kmer_DF$TEclass  <- TE_GTF_defrag$class_id[match(rownames(Kmer_DF),TE_GTF_defrag$defrag_transcript_id)]
	Kmer_DF$TEfam    <- TE_GTF_defrag$family_id[match(rownames(Kmer_DF),TE_GTF_defrag$defrag_transcript_id)]
	Kmer_DF$TEsubfam <- TE_GTF_defrag$gene_id[match(rownames(Kmer_DF),TE_GTF_defrag$defrag_transcript_id)]
	Kmer_DF$TEloci   <- TE_GTF_defrag$defrag_transcript_id[match(rownames(Kmer_DF),TE_GTF_defrag$defrag_transcript_id)]
	#
	if(TE_type=="TEclass"){
		TE_list <- unique(Kmer_DF$TEclass)
	} else if(TE_type=="TEfam"){
		TE_list <- unique(Kmer_DF$TEfam)
	} else if(TE_type=="TEsubfam"){
		TE_list <- unique(Kmer_DF$TEsubfam)
	} else{
		stop("TE_type argument must be \"TEclass\", \"TEfam\" or \"TEsubfam\".")
	}
	# Use get_contTable_vals find values for each TE
	DF_TE_fisher <- data.frame(TE=rep(TE_list, each=N_clusters),
								Cluster=rep( LETTERS[1:N_clusters], length(TE_list)) )
	DF_TE_fisher <- as.data.frame( t( apply(DF_TE_fisher,1, function(x) c(x, get_contTable_vals(TE_type, x[1], x[2], Kmer_DF) )) ) )
	# Run fisher.test for each line with apply
	DF_TE_fisher$P.value   <- apply( DF_TE_fisher,1, function(x) fisher.test(matrix(as.numeric(x[3:6]),nr=2))$p.value )
	DF_TE_fisher$OddRatio  <- apply( DF_TE_fisher,1, function(x) fisher.test(matrix(as.numeric(x[3:6]),nr=2))$estimate )
	DF_TE_fisher$Padjust   <- p.adjust(DF_TE_fisher$P.value, method="fdr")
	DF_TE_fisher$Log10_FDR <- log10(DF_TE_fisher$Padjust)*-1
	DF_TE_fisher$Log2_OR   <- log2(DF_TE_fisher$OddRatio)
	#
	return(DF_TE_fisher)
}

### Making unique TEclass/TEfamily/TEsubfam table
tic("Making unique TEclass/TEfamily/TEsubfam table")
TE_table <- unique(data.frame(gene_id=TE_GTF_defrag$gene_id, family_id=TE_GTF_defrag$family_id, class_id=TE_GTF_defrag$class_id))
toc()

### Run enrichment
tic("Run enrichment")
TEclass_fisher   <- fisher_test_cluster_TEloci(Kmean_clust_rename, 7, "TEclass")
TEfam_fisher     <- fisher_test_cluster_TEloci(Kmean_clust_rename, 7, "TEfam")
TEfam_sig        <- TEfam_fisher$TE[which(TEfam_fisher$Padjust<0.05)]
TEfam_fisher_sub <- TEfam_fisher[which(TEfam_fisher$TE %in% TEfam_sig),]
# Join TE class and TE family
TEclass_fisher$TE <- paste0("C:",TEclass_fisher$TE)
TEclass_fam_fisher <- rbind(TEclass_fisher, TEfam_fisher_sub)
toc()

### Plot
tic("Plotting")
# Change Inf and -Inf from significant lines with the max and min log2 OR value
MAX_Log2_OR <- max(TEclass_fam_fisher$Log2_OR[is.finite(TEclass_fam_fisher$Log2_OR)])
MIN_Log2_OR <- min(TEclass_fam_fisher$Log2_OR[is.finite(TEclass_fam_fisher$Log2_OR)])
#
TEclass_fam_fisher$Log2_OR[which(TEclass_fam_fisher$Log2_OR==Inf & TEclass_fam_fisher$Padjust<0.05)]  <- MAX_Log2_OR
TEclass_fam_fisher$Log2_OR[which(TEclass_fam_fisher$Log2_OR==-Inf & TEclass_fam_fisher$Padjust<0.05)] <- MIN_Log2_OR
# p-value breaks
fdr_breaks <- c(0, log10(0.05)*-1, 20, 50, 100 )
fdr_labels <- c("Padj. > 0.05","Padj. < 0.05", paste0("Padj. < 1e-",fdr_breaks[3:length(fdr_breaks)]))
max_log10_FDR <- max(round(TEclass_fam_fisher$Log10_FDR[which(abs(TEclass_fam_fisher$Log10_FDR)!=Inf)]))
for (fdr in fdr_breaks) {
	if(fdr == fdr_breaks[1]) {
		TEclass_fam_fisher$Log10_FDR_breaks <- max_log10_FDR
		TEclass_fam_fisher$Log10_FDR_breaks[which(TEclass_fam_fisher$Log10_FDR<fdr_breaks[2])] <- fdr
	} else if (fdr == fdr_breaks[length(fdr_breaks)]) {
		TEclass_fam_fisher$Log10_FDR_breaks[which(TEclass_fam_fisher$Log10_FDR>=fdr)] <- fdr
	} else{
		fdr_idx <- which(fdr_breaks==fdr)
		TEclass_fam_fisher$Log10_FDR_breaks[which(TEclass_fam_fisher$Log10_FDR<fdr_breaks[fdr_idx+1] & TEclass_fam_fisher$Log10_FDR>=fdr )] <- fdr
	}
}
# Order clusters
TEclass_fam_fisher$Cluster <- factor(TEclass_fam_fisher$Cluster, levels=rev(unique(TEclass_fam_fisher$Cluster)))
# Order TE classes and TE families
sorted_classes <- sort(unique(TE_table$class_id))
sorted_fams <- c()
for (C in sorted_classes) {
	sorted_fams <- c(sorted_fams, unique(sort(as.character(TE_table$family_id[which(TE_table$class_id==C)]))) )
}
sorted_fams_sub <- sorted_fams[which(sorted_fams %in% unique(as.character(TEclass_fam_fisher$TE)))]
sorted_fams_sub_append_class <- c(paste0("C:",sorted_classes), sorted_fams_sub)
TEclass_fam_fisher$TE <- factor(TEclass_fam_fisher$TE, levels=sorted_fams_sub_append_class )
# Create a new column of significant vs non-significant
TEclass_fam_fisher$Padjust_sig <- rep("Not_significant",nrow(TEclass_fam_fisher))
TEclass_fam_fisher$Padjust_sig[which(TEclass_fam_fisher$Padjust<0.05)] <- "Significant"
# TE class and family separate
TEclass_fam_fisher_subC <- TEclass_fam_fisher[grep("C:",TEclass_fam_fisher$TE),]
TEclass_fam_fisher_subF <- TEclass_fam_fisher[grep("C:",TEclass_fam_fisher$TE, invert=TRUE),]
TEclass_fam_fisher_subC_lev <- levels(TEclass_fam_fisher_subC$TE)[grep("C:",levels(TEclass_fam_fisher_subC$TE))]
TEclass_fam_fisher_subF_lev <- levels(TEclass_fam_fisher_subF$TE)[grep("C:",levels(TEclass_fam_fisher_subF$TE), invert=TRUE)]
TEclass_fam_fisher_subC$TE <- factor(as.character(TEclass_fam_fisher_subC$TE), levels=TEclass_fam_fisher_subC_lev)
TEclass_fam_fisher_subF$TE <- factor(as.character(TEclass_fam_fisher_subF$TE), levels=TEclass_fam_fisher_subF_lev)
#
TE_dotplot_subC <- ggplot(TEclass_fam_fisher_subC, aes(TE, Cluster)) +
    geom_point(aes(color=Log2_OR, size=Log10_FDR_breaks)) +
    scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", limits=c(MIN_Log2_OR, MAX_Log2_OR)) +
    ylab("Clusters") +
    theme(plot.title=element_text(hjust=0.5)) + theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position="none") + 
    scale_size_continuous(breaks=fdr_breaks, range=c(1,7), labels=fdr_labels)
TE_dotplot_subF <- ggplot(TEclass_fam_fisher_subF, aes(TE, Cluster)) +
    geom_point(aes(color=Log2_OR, size=Log10_FDR_breaks)) +
    scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", limits=c(MIN_Log2_OR, MAX_Log2_OR)) + 
    ylab("Clusters") +
    theme(plot.title=element_text(hjust=0.5)) + theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank()) +
    scale_size_continuous(breaks=fdr_breaks, range=c(1,7), labels=fdr_labels)
# Put one next to the other
library(patchwork)
TE_dotplot_sub_join <- TE_dotplot_subC + TE_dotplot_subF + plot_layout(ncol=2, widths=c(1,5))
# Save
OUT_DIR <- "./data/TE_enrichment"
pdf(file.path(OUT_DIR,"TE_enrich_TEclass_n_fam_min10.pdf"), height=5,width=10)
print(TE_dotplot_sub_join)
dev.off()
toc()

############################