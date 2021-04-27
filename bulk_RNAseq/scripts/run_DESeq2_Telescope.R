#!/usr/bin/env Rscript
#run_DESeq2_Telescope.R
# Perform DE analysis of TE loci defrag with 
# the count data from telescope stranded.

suppressMessages(library(tictoc))
suppressMessages(library(DESeq2))
suppressMessages(library(reshape2))
suppressMessages(library(gridExtra))
suppressMessages(library(dplyr))
suppressMessages(library(RUVSeq))


## Read telescope counts
TELES_DIR <- "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/merge_Teles_counts"
writeLines("Reading telescope stranded results")
tic("  - Reading telescope stranded results")
Telescope_count_fwd_fwd <- read.delim("./data/telescope/fwd_fwd/Edited_telescope_report/join_telesRep.final_count.dfrag.tab", header=TRUE)
Telescope_count_rev_rev <- read.delim("./data/telescope/rev_rev/Edited_telescope_report/join_telesRep.final_count.dfrag.tab", header=TRUE)
rownames(Telescope_count_fwd_fwd) <- Telescope_count_fwd_fwd$TEs; Telescope_count_fwd_fwd <- Telescope_count_fwd_fwd[,2:ncol(Telescope_count_fwd_fwd)]
rownames(Telescope_count_rev_rev) <- Telescope_count_rev_rev$TEs; Telescope_count_rev_rev <- Telescope_count_rev_rev[,2:ncol(Telescope_count_rev_rev)]
colnames(Telescope_count_fwd_fwd) <- gsub("^X","",colnames(Telescope_count_fwd_fwd))
colnames(Telescope_count_rev_rev) <- gsub("^X","",colnames(Telescope_count_rev_rev))
toc()
#Merge FWD_FWD and REV_REV
writeLines("Merging forward and reverse data")
tic("  - Merging forward and reverse data")
Telescope_count <- Telescope_count_fwd_fwd+Telescope_count_rev_rev
toc()

##Filter: filter out non-expressed TEs, by requiring more 
#        than 5 reads in at least two samples for each gene.
writeLines("Filtering out non-expressed TEs (>5 reads in at least 2 samples)")
tic("  - Filtering")
filter <- apply(Telescope_count, 1, function(x) length(x[x>5])>=2)
Telescope_count <- Telescope_count[filter,]
toc()

## Join TE counts and gene counts (use TEtranscripts gene counts)
#   This is necessary because genes contribute more to the variance in
#   the samples. Therefore, running DESeq() without gene counts might 
#   lead into not accurate estimation of size factors and dispersion.
writeLines("Join TE counts and gene counts for a better size factors and dispersion estimation")
tic("  - Read dds_TEtrans_GE.rds"); dds_TEtrans_GE <- readRDS("./data/TEtranscripts/dds_TEtrans_GE.rds"); toc()
if (all(colnames(dds_TEtrans_GE)==colnames(Telescope_count)) == FALSE){
	stop("ERROR: column names are different between dds_TEtrans_GE and TE counts")
}
Telescope_count <- rbind(counts(dds_TEtrans_GE,normalize=FALSE),
                                Telescope_count)


## Calculate RUV factor normalization
# Read TEtranscripts counts
counts_path <- "./data/TEtranscripts/join.cntTable"
tic("Reading counts multi")
counts <- read.delim(counts_path, header=TRUE, row.names="Feature")
toc()
colnames(counts) <- gsub("^X","",colnames(counts))
# Re-order replicate columns
stages <- unique(gsub("_rep_[1-5]","",colnames(counts)))
col_order <- c()
for (s in stages) {
        idx <- order(colnames(counts)[grep(s, colnames(counts))])
        col_order <- c(col_order, colnames(counts)[grep(s, colnames(counts))][idx])
}
counts <- counts[match(col_order, colnames(counts))]
#Filter like the TE counts
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
filtered <- counts[filter,]
# Separate genes and ercc lines
genes <- filtered[grep("^ENSDARG", rownames(filtered)),]
ercc  <- filtered[grep("^ERCC", rownames(filtered)),]
#
x <- factor(gsub("_rep_[1-5]","",colnames(genes)), levels=unique(gsub("_rep_[1-5]","",colnames(genes))) )
set <- newSeqExpressionSet(as.matrix(filtered), phenoData=data.frame(x, row.names=colnames(filtered)))
# Calculate the factors of unwanted variation
set1 <- RUVg(set, rownames(ercc), k=1)
RUV_factor <- pData(set1)
RUV_factor_DF <- data.frame(Sample_name=RUV_factor$x, Unwanted_varaition=RUV_factor$W_1)
RUV_factor_DF$Samples <- factor(rownames(RUV_factor), levels=as.character(rownames(RUV_factor)))


## Load counts to DESeq2
# Make a coldata file
coldata <- data.frame(Sample_name=colnames(Telescope_count),
                      Stage=gsub("_rep_[1-5]","",colnames(Telescope_count)),
                      Replicate=gsub(".*_rep_","Rep_",colnames(Telescope_count)) )
coldata$RUV_factor <- RUV_factor_DF$Unwanted_varaition
# Load
writeLines("Loading counts to DESeq2")
tic("  - Loading counts to DESeq2")
dds_TEloci_defrag <- DESeqDataSetFromMatrix(countData = Telescope_count,
                                            colData = coldata,
                                            design= ~ RUV_factor + Stage)
toc()

## Run DESeq for gene counts
writeLines("Running DESeq2")
tic("  - Running DESeq2")
dds_TEloci_defrag <- DESeq(dds_TEloci_defrag)
toc()

## Discard gene lines from the dds object
writeLines("Discard gene lines from DESeq objects")
dds_TEloci_defrag <- dds_TEloci_defrag[grep("^ENSDARG", rownames(dds_TEloci_defrag), invert=TRUE),]

## Save dds
writeLines("Saving ./data/DE_analysis/dds_TEloci_defrag.rds")
tic("  - Saving dds_TEloci_defrag.rds")
saveRDS(dds_TEloci_defrag, "./data/DE_analysis/dds_TEloci_defrag.rds")
toc()
## Reads dds
#tic("Read dds_TEloci_defrag.rds")
#dds_TEloci_defrag <- readRDS(file.path(DE_TE_DIR,"./data/DE_analysis/dds_TEloci_defrag.rds"))
#toc()

###############