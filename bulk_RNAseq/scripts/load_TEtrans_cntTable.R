#!/usr/bin/env Rscript
#load_TEtrans_cntTable.R

suppressMessages(library(tictoc))
suppressMessages(library(DESeq2))

## Read TEtranscripts counts
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

## Remove ERCC lines
counts <- counts[grep("^ERCC", rownames(counts), invert=TRUE),]

## Load counts to DESeq2
coldata <- data.frame(Sample_name=colnames(counts),
                      Stage=gsub("_rep_[1-5]","",colnames(counts)),
                      Replicate=gsub(".*_rep_","Rep_",colnames(counts)) )
# Load
tic("Load counts to DESeq2")
dds_TEtrans <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = coldata,
                                      design = ~Stage)
toc()

## Separate genes and TEs
dds_TEtrans_GE <- dds_TEtrans[grep("^ENSDARG",rownames(dds_TEtrans)),]
dds_TEtrans_TE <- dds_TEtrans[grep("^ENSDARG",rownames(dds_TEtrans), invert=TRUE),]

## Save RDS
tic("Saving dds_TEtrans_GE.rds")
saveRDS(dds_TEtrans_GE, "./data/TEtranscripts/dds_TEtrans_GE.rds")
toc()
tic("Saving dds_TEtrans_TE.rds")
saveRDS(dds_TEtrans_TE, "./data/TEtranscripts/dds_TEtrans_TE.rds")
toc()
#
writeLines("*** Done ***")

#####################