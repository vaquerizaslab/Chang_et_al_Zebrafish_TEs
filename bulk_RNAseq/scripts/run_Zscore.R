#!/usr/bin/env Rscript
#run_Zscore.R
# Load Telescope counts.
# Extract normalized count matrix of the DE TE loci.
# Transform matrix to Z-score.

suppressMessages(library(tictoc))
tic("Loading packages")
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
toc()

### Load counts for genes and TEloci
tic("Read dds_TEloci_defrag.rds")
dds_TEloci_defrag <- readRDS("./data/DE_analysis/dds_TEloci_defrag.rds")
toc()

### Load TE annotation with read-through categories
### Load TE annotations with extended 3' UTR data
tic("Loading TE annotation")
TE_GTF_DFRAG_path <- "./data/danRer11.TEtrans_uID.dfrag.classified.gtf"
TE_GTF_defrag <- import(TE_GTF_DFRAG_path)
toc(); writeLines("")

### Subset non read-through TE loci
TE_GTF_defrag_selfExp <- TE_GTF_defrag[which(TE_GTF_defrag$Rt_catg=="Self_expressed")]

### Extract normalized counts
counts_TEl_defrag <- counts(dds_TEloci_defrag, normalized=TRUE)

### Reduce to DE TEs
sig_DE_logFC0_padj01 <- read.table("./data/DE_analysis/pairwise_results/sig_DE/All_sigDE_TEloci.pAdj0.01.log2FC0.txt", header=FALSE)
sig_DE_logFC0_padj01 <- sig_DE_logFC0_padj01$V1
counts_TEl_defrag_DE_all <- counts_TEl_defrag[which(rownames(counts_TEl_defrag) %in% sig_DE_logFC0_padj01),]
counts_TEl_defrag_DE_selfExp <- counts_TEl_defrag_DE_all[which(rownames(counts_TEl_defrag_DE_all) %in% TE_GTF_defrag_selfExp$defrag_transcript_id),]

### Remove TEs with less than 10 reads in any stage
counts_TEl_defrag_DE_selfExp_min10 <- counts_TEl_defrag_DE_selfExp[apply(counts_TEl_defrag_DE_selfExp, 1, function(x) !(all(x<10))),]

### Transform to Z-score
tic("Transform to Z-score")
Zscore_TEl_defrag_DE_selfExp_min10 <- t( apply(counts_TEl_defrag_DE_selfExp_min10, 1, function(x) (x-mean(x))/sd(x) ) )
toc()

### Save Zscore matrix
tic("Save Z-score matrices")
OUT_DIR <- "./data/Zscore_matrix"
saveRDS(Zscore_TEl_defrag_DE_selfExp_min10, file.path(OUT_DIR,"Zscore_TEl_defrag_DE_selfExp.min10.rds"))
toc()

################