#!/usr/bin/env Rscript
# get_extended_3UTR.R
# R script to extract extended 3' UTR
# regions in novel isoforms of known genes.

suppressWarnings(suppressMessages(library(optparse)))
option_list <- list(
  make_option(c("-s","--stringtie"), type="character", default=NULL, 
              help="Path to stringtie GTF output."),
  make_option(c("-r","--reference"), type="character", default=NULL, 
              help="Path to reference GTF."),
  make_option(c("-o","--output"), type="character", default=NULL, 
              help="Path BED output.")
);
opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

#### Print parameters
writeLines( paste0("Parameter --stringtie: ", opt$stringtie) )
writeLines( paste0("Parameter --reference: ", opt$reference) )
writeLines( paste0("Parameter --output: ", opt$output) )

### Load packages
suppressWarnings(suppressMessages(library(tictoc)))
tic("Loading packages")
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(regioneR)))
toc()

### Read GTF files
tic("Reading GTFs")
STR_GTF <- import(opt$stringtie)
GE_GTF <- import(opt$reference)
toc()

### Extract exons
tic("Extracting exons and transcript")
STR_GTF_exons <- STR_GTF[which(STR_GTF$type=="exon")]
GE_GTF_exons <- GE_GTF[which(GE_GTF$type=="exon")]
GE_GTF_trans <- GE_GTF[which(GE_GTF$type=="transcript")]
toc()

### Select last exon for each gene
tic("Selection last exon for each gene")
GE_GTF_exons_fwd <- GE_GTF_exons[which(strand(GE_GTF_exons)=="+")]
GE_GTF_exons_rev <- GE_GTF_exons[which(strand(GE_GTF_exons)=="-")]
GE_GTF_last_exon_fwd <- as.data.frame(GE_GTF_exons_fwd) %>% group_by(gene_id) %>% slice(which.max(end))
GE_GTF_last_exon_rev <- as.data.frame(GE_GTF_exons_rev) %>% group_by(gene_id) %>% slice(which.min(start))
GE_GTF_last_exon_fwd <- makeGRangesFromDataFrame(GE_GTF_last_exon_fwd, keep.extra.columns=TRUE)
GE_GTF_last_exon_rev <- makeGRangesFromDataFrame(GE_GTF_last_exon_rev, keep.extra.columns=TRUE)
GE_GTF_last_exon_fwd <- sort(GE_GTF_last_exon_fwd)
GE_GTF_last_exon_rev <- sort(GE_GTF_last_exon_rev)
toc()

### Subset stringtie exons that overlap last exons
tic("Subset Stringtie exons that overlap last exons")
STR_GTF_exons_overlap_fwd <- STR_GTF_exons[unique(queryHits(findOverlaps(STR_GTF_exons, GE_GTF_last_exon_fwd)))]
STR_GTF_exons_overlap_rev <- STR_GTF_exons[unique(queryHits(findOverlaps(STR_GTF_exons, GE_GTF_last_exon_rev)))]
toc()

### Cut only the fraction that is downstream the reference 3' end.
tic("Extract only extended 3' UTR transcription regions")
GE_GTF_trans_fwd <- GE_GTF_trans[which(strand(GE_GTF_trans)=="+")]
GE_GTF_trans_rev <- GE_GTF_trans[which(strand(GE_GTF_trans)=="-")]
extended_regions_fwd <- subtractRegions(STR_GTF_exons_overlap_fwd, GE_GTF_trans_fwd)
extended_regions_rev <- subtractRegions(STR_GTF_exons_overlap_rev, GE_GTF_trans_rev)
toc()

### Add name of the closest upstream gene for each pervasive region
tic("Pair extended transcription with its nearest gene")
extended_regions_fwd$gene_id      <- GE_GTF_last_exon_fwd$gene_id[nearest(extended_regions_fwd, GE_GTF_last_exon_fwd)]
extended_regions_fwd$gene_name    <- GE_GTF_last_exon_fwd$gene_name[nearest(extended_regions_fwd, GE_GTF_last_exon_fwd)]
extended_regions_fwd$gene_biotype <- GE_GTF_last_exon_fwd$gene_biotype[nearest(extended_regions_fwd, GE_GTF_last_exon_fwd)]
extended_regions_rev$gene_id      <- GE_GTF_last_exon_rev$gene_id[nearest(extended_regions_rev, GE_GTF_last_exon_rev)]
extended_regions_rev$gene_name    <- GE_GTF_last_exon_rev$gene_name[nearest(extended_regions_rev, GE_GTF_last_exon_rev)]
extended_regions_rev$gene_biotype <- GE_GTF_last_exon_rev$gene_biotype[nearest(extended_regions_rev, GE_GTF_last_exon_rev)]
toc()

### Merge + and -
extended_regions <- c(extended_regions_fwd, extended_regions_rev)
extended_regions <- sort(extended_regions, ignore.strand=TRUE)

### Save in BED format
tic("Save extended 3' UTR transcription as BED file")
extended_regions_BED <- data.frame(seqnames=seqnames(extended_regions),
  starts=start(extended_regions)-1,
  ends=end(extended_regions),
  gene_id=extended_regions$gene_id,
  scores=c(rep(".", length(extended_regions))),
  strands=strand(extended_regions),
  gene_name=extended_regions$gene_name,
  gene_biotype=extended_regions$gene_biotype)
write.table(extended_regions_BED, file=opt$output, quote=F, sep="\t", row.names=F, col.names=F)
toc()

##############