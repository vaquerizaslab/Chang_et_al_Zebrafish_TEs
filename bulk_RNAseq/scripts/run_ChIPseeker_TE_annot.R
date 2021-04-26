#!/usr/bin/env Rscript
#run_ChIPseeker_TE_annot.R
# Use ChIPseeker to annotate the position of TE with respect to gene annotations

suppressWarnings(suppressMessages(library(optparse)))
option_list <- list(
  make_option(c("-t","--te_gtf"), type="character", default=NULL, 
              help="Path to TE GTF."),
  make_option(c("-g","--gene_gtf"), type="character", default=NULL, 
              help="Path to gene GTF (protein coding genes subset)."),
  make_option(c("-o","--output"), type="character", default=NULL, 
              help="Path to output GTF file.")
);
opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

### #Load packages
suppressMessages(library(tictoc))
tic("Loading packages")
suppressMessages(library(DESeq2))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
toc(); writeLines("")

### Load TE annotations with defrag TEs
tic("Loading TE defrag GTF")
TE_GTF_defrag <- import(opt$te_gtf)
toc(); writeLines("")

### Build TxDb for GRCz11 Protein coding genes
tic("Making TxDb for GRCz11 (PC genes)")
GE_GTF_PC_path <- opt$gene_gtf
chrom_size_p <- "" #Add here the path to a 2 column file with chr in the 1st column and the number of bp in the 2nd column
chrom_size <- read.delim(chrom_size_p, header=FALSE, col.names=c("chrom", "length"))
chrom_size$is_circular <- rep(FALSE, nrow(chrom_size))
TxDb_GRCz11_PC <- makeTxDbFromGFF(file=GE_GTF_PC_path,
                                   chrominfo=chrom_size,
                                   dataSource="ensemblgenomes",
                                   organism="Danio rerio")
toc()

### Save TxDb for GRCz11 (PC genes)
tic("Save TxDb for GRCz11 (PC genes)")
TxDb_GRCz11_PC_path <- sub(".gtf",".sqlite",opt$gene_gtf)
saveDb(TxDb_GRCz11_PC, TxDb_GRCz11_PC_path)
toc()

### Annotate with ChIPseeker 
tic("ChIPseeker annotate all TEs vs PC genes")
peakAnno_TE_PCg <- annotatePeak(TE_GTF_defrag, tssRegion=c(-3000, 3000), TxDb=TxDb_GRCz11_PC)
peakAnno_TE_PCg_GR <- as.GRanges(peakAnno_TE_PCg)
peakAnno_TE_PCg_order <- annotatePeak(TE_GTF_defrag, tssRegion=c(-3000, 3000), TxDb=TxDb_GRCz11_PC, genomicAnnotationPriority=c("5UTR","3UTR","Exon","Promoter","Intron","Downstream","Intergenic"))
peakAnno_TE_PCg_order_GR <- as.GRanges(peakAnno_TE_PCg_order)
# SS = Same Strand
# OS = Opposite Strand
#Fix missing Promoter annotations
intersect_Inter_Prom <- intersect(grep("Distal Intergenic",peakAnno_TE_PCg_order_GR$annotation), grep("Promoter", peakAnno_TE_PCg_GR$annotation))
intersect_Downs_Prom <- intersect(grep("Downstream",peakAnno_TE_PCg_order_GR$annotation), grep("Promoter", peakAnno_TE_PCg_GR$annotation))
peakAnno_TE_PCg_order_GR$annotation_fix <- peakAnno_TE_PCg_order_GR$annotation
peakAnno_TE_PCg_order_GR$annotation_fix[intersect_Inter_Prom] <- peakAnno_TE_PCg_GR$annotation[intersect_Inter_Prom]
peakAnno_TE_PCg_order_GR$annotation_fix[intersect_Downs_Prom] <- peakAnno_TE_PCg_GR$annotation[intersect_Downs_Prom]
# Fix v3: If distanceToTSS is >0, then it shouldn't be a Promoter
prom_distTSS_pos_idx <- which((peakAnno_TE_PCg_order_GR$annotation_fix=="Promoter (<=1kb)" | peakAnno_TE_PCg_order_GR$annotation_fix=="Promoter (1-2kb)" | peakAnno_TE_PCg_order_GR$annotation_fix=="Promoter (2-3kb)") & peakAnno_TE_PCg_order_GR$distanceToTSS>0)
peakAnno_TE_PCg_order_GR$annotation_fix[prom_distTSS_pos_idx] <- peakAnno_TE_PCg_order_GR$annotation[prom_distTSS_pos_idx]
#
peakAnno_TE_PCg_order_GR$TE_GE_strand <- rep("SS", length(peakAnno_TE_PCg_order_GR))
peakAnno_TE_PCg_order_GR$TE_GE_strand[which(peakAnno_TE_PCg_order_GR$geneStrand==2 & strand(peakAnno_TE_PCg_order_GR)=="+")] <- "OS"
peakAnno_TE_PCg_order_GR$TE_GE_strand[which(peakAnno_TE_PCg_order_GR$geneStrand==1 & strand(peakAnno_TE_PCg_order_GR)=="-")] <- "OS"
peakAnno_TE_PCg_order_GR$annotation_2 <- gsub("Exon .*","Exon",gsub("Intron .*","Intron",peakAnno_TE_PCg_order_GR$annotation_fix))
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="3' UTR")] <- "3_UTR"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="5' UTR")] <- "5_UTR"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Distal Intergenic")] <- "Intergenic"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Downstream (<1kb)")] <- "Downstream_1kb"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Downstream (1-2kb)")] <- "Downstream_2kb"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Downstream (2-3kb)")] <- "Downstream_3kb"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Promoter (<=1kb)")] <- "Promoter_1kb"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Promoter (1-2kb)")] <- "Promoter_2kb"
peakAnno_TE_PCg_order_GR$annotation_2[which(peakAnno_TE_PCg_order_GR$annotation_2=="Promoter (2-3kb)")] <- "Promoter_3kb"
peakAnno_TE_PCg_order_GR$annotation_2 <- paste0(peakAnno_TE_PCg_order_GR$annotation_2,".",peakAnno_TE_PCg_order_GR$TE_GE_strand)
peakAnno_TE_PCg_order_GR$annotation_3 <- gsub("_[123]kb","",peakAnno_TE_PCg_order_GR$annotation_2)
# Subset important info (include transcript that the annotation feature refers)
peakAnno_TE_PCg_GR_sub <- peakAnno_TE_PCg_order_GR[,c(1:8,15:16,18:21)]
toc()

### Save as GTF
tic(paste0("Exporting: ",opt$output))
export(peakAnno_TE_PCg_GR_sub, opt$output))
toc()


##############