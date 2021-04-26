#!/usr/bin/env Rscript
# combine_ext_3UTRs.R
# Load extended 3' UTRs from each stage and combine it into a single set of regions.
# Tracks made by bamCoverage will be used to discard regions that have 0 signal on them (potential,
# false positives from stringtie).

suppressWarnings(suppressMessages(library(optparse)))
option_list <- list(
  make_option(c("-t","--te_gtf"), type="character", default=NULL, 
              help="Path to TE GTF."),
  make_option(c("-d","--dir_ext_UTR"), type="character", default=NULL, 
              help="Path to directory with extended 3' UTR files (they must have suffix like this *.ext_3UTR.CPM.bed)."),
  make_option(c("-o","--output"), type="character", default=NULL, 
              help="Path to output GTF file.")
);
opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

### Load packages
suppressMessages(library(tictoc))
tic("Loading packages")
suppressMessages(library(rtracklayer))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(UpSetR))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
toc(); writeLines("")


### Load TE annotations
tic("Loading TE annotations")
TE_GTF <- import(opt$te_gtf)
TE_GTF_sub <- TE_GTF[,c(3,5:8,14)]
toc()

### Load pervasive transcription (PT) regions
DIR_EXT_3UTR <- opt$dir_ext_UTR
all_ext_3UTR_paths <- c(file.path(DIR_EXT_3UTR,"1_cell.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"2_cell.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"128_cell.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"1k_cell.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"dome.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"50pc_epiboly.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"shield.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"75pc_epiboly.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"1_4_somites.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"14_19_somites.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"20_25_somites.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"prim_5.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"prim_15.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"prim_25.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"long_pec.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"protruding_mouth.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"day_4.ext_3UTR.CPM.bed"),
                        file.path(DIR_EXT_3UTR,"day_5.ext_3UTR.CPM.bed"))
tic("Loading pervasive transcription regions")
EXT_1_cell           <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[1],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_2_cell           <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[2],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_128_cell         <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[3],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_1k_cell          <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[4],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_dome             <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[5],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_50pc_epiboly     <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[6],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_shield           <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[7],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_75pc_epiboly     <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[8],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_1_4_somites      <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[9],  header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_14_19_somites    <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[10], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_20_25_somites    <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[11], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_prim_5           <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[12], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_prim_15          <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[13], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_prim_25          <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[14], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_long_pec         <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[15], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_protruding_mouth <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[16], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_day_4            <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[17], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
EXT_day_5            <- makeGRangesFromDataFrame(read.delim( all_ext_3UTR_paths[18], header=FALSE, col.names=c("chr","start","end","gene_id","score","strand","gene_name","gene_biotype","CPM_mean","CPM_max") ), keep.extra.columns=TRUE)
#
EXT_stages <- GRangesList(`1_cell`=EXT_1_cell,
                          `2_cell`=EXT_2_cell,
                          `128_cell`=EXT_128_cell,
                          `1k_cell`=EXT_1k_cell,
                          `dome`=EXT_dome,
                          `50pc_epiboly`=EXT_50pc_epiboly,
                          `shield`=EXT_shield,
                          `75pc_epiboly`=EXT_75pc_epiboly,
                          `1_4_somites`=EXT_1_4_somites,
                          `14_19_somites`=EXT_14_19_somites,
                          `20_25_somites`=EXT_20_25_somites,
                          `prim_5`=EXT_prim_5,
                          `prim_15`=EXT_prim_15,
                          `prim_25`=EXT_prim_25,
                          `long_pec`=EXT_long_pec,
                          `protruding_mouth`=EXT_protruding_mouth,
                          `day_4`=EXT_day_4,
                          `day_5`=EXT_day_5)
stages <- c("1_cell","2_cell","128_cell","1k_cell","dome","50pc_epiboly","shield","75pc_epiboly","1_4_somites","14_19_somites","20_25_somites","prim_5","prim_15","prim_25","long_pec","protruding_mouth","day_4","day_5")
# Remove unnecessary objects
rm(EXT_1_cell); rm(EXT_2_cell); rm(EXT_128_cell); rm(EXT_1k_cell); rm(EXT_dome); rm(EXT_50pc_epiboly); rm(EXT_shield); rm(EXT_75pc_epiboly); rm(EXT_1_4_somites); rm(EXT_14_19_somites); rm(EXT_20_25_somites); rm(EXT_prim_5); rm(EXT_prim_15); rm(EXT_prim_25); rm(EXT_long_pec); rm(EXT_protruding_mouth); rm(EXT_day_4); rm(EXT_day_5)
# Discard MT perv trans regions
EXT_stages_noMT <- GRangesList(`1_cell`=EXT_stages$`1_cell`[which(seqnames(EXT_stages$`1_cell`)!="chrM")],
                               `2_cell`=EXT_stages$`2_cell`[which(seqnames(EXT_stages$`2_cell`)!="chrM")],
                               `128_cell`=EXT_stages$`128_cell`[which(seqnames(EXT_stages$`128_cell`)!="chrM")],
                               `1k_cell`=EXT_stages$`1k_cell`[which(seqnames(EXT_stages$`1k_cell`)!="chrM")],
                               `dome`=EXT_stages$`dome`[which(seqnames(EXT_stages$`dome`)!="chrM")],
                               `50pc_epiboly`=EXT_stages$`50pc_epiboly`[which(seqnames(EXT_stages$`50pc_epiboly`)!="chrM")],
                               `shield`=EXT_stages$`shield`[which(seqnames(EXT_stages$`shield`)!="chrM")],
                               `75pc_epiboly`=EXT_stages$`75pc_epiboly`[which(seqnames(EXT_stages$`75pc_epiboly`)!="chrM")],
                               `1_4_somites`=EXT_stages$`1_4_somites`[which(seqnames(EXT_stages$`1_4_somites`)!="chrM")],
                               `14_19_somites`=EXT_stages$`14_19_somites`[which(seqnames(EXT_stages$`14_19_somites`)!="chrM")],
                               `20_25_somites`=EXT_stages$`20_25_somites`[which(seqnames(EXT_stages$`20_25_somites`)!="chrM")],
                               `prim_5`=EXT_stages$`prim_5`[which(seqnames(EXT_stages$`prim_5`)!="chrM")],
                               `prim_15`=EXT_stages$`prim_15`[which(seqnames(EXT_stages$`prim_15`)!="chrM")],
                               `prim_25`=EXT_stages$`prim_25`[which(seqnames(EXT_stages$`prim_25`)!="chrM")],
                               `long_pec`=EXT_stages$`long_pec`[which(seqnames(EXT_stages$`long_pec`)!="chrM")],
                               `protruding_mouth`=EXT_stages$`protruding_mouth`[which(seqnames(EXT_stages$`protruding_mouth`)!="chrM")],
                               `day_4`=EXT_stages$`day_4`[which(seqnames(EXT_stages$`day_4`)!="chrM")],
                               `day_5`=EXT_stages$`day_5`[which(seqnames(EXT_stages$`day_5`)!="chrM")])
toc(); writeLines("")

### PT filter regions with mean CPM == 0
CPM_bt_0_idx  <- lapply(EXT_stages_noMT, function(x) x$CPM_mean > 0)
sub_GRlist <- function(GR_list, sub_list){
    #Check that both lists have the same length
    if( length(GR_list) != length(sub_list) ){
        writeLines("ERROR: GR list and subset list must have the same length.")
        break
    }
    #Check that the lengths of the lists inside GR_list and sub_list are the same
    if( all(unlist(lapply(GR_list, length)) == unlist(lapply(sub_list, length)))==FALSE ){
        writeLines("ERROR: GR list and subset list vectors must have the same length.")
        break
    }
    #Subset
    GR_list_sub <- GR_list
    for ( i in seq(1,length(GR_list)) ) {
        GR_list_sub[[i]] <- GR_list[[i]][sub_list[[i]]]
    }
    #Return
    return(GR_list_sub)
}
#
EXT_stages_noMT_CPM_0  <- sub_GRlist(EXT_stages_noMT, CPM_bt_0_idx)

### Get TEs overlapping PT
TE_GTF_sub_EXT <- TE_GTF_sub
TE_GTF_sub_EXT$EXT_1_cell           <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_2_cell           <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_128_cell         <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_1k_cell          <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_dome             <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_50pc_epiboly     <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_shield           <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_75pc_epiboly     <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_1_4_somites      <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_14_19_somites    <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_20_25_somites    <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_prim_5           <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_prim_15          <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_prim_25          <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_long_pec         <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_protruding_mouth <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_day_4            <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_day_5            <- rep("No",length(TE_GTF_sub_EXT))
TE_GTF_sub_EXT$EXT_1_cell[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="1_cell")])) ]                     <- "Yes"
TE_GTF_sub_EXT$EXT_2_cell[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="2_cell")])) ]                     <- "Yes"
TE_GTF_sub_EXT$EXT_128_cell[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="128_cell")])) ]                 <- "Yes"
TE_GTF_sub_EXT$EXT_1k_cell[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="1k_cell")])) ]                   <- "Yes"
TE_GTF_sub_EXT$EXT_dome[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="dome")])) ]                         <- "Yes"
TE_GTF_sub_EXT$EXT_50pc_epiboly[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="50pc_epiboly")])) ]         <- "Yes"
TE_GTF_sub_EXT$EXT_shield[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="shield")])) ]                     <- "Yes"
TE_GTF_sub_EXT$EXT_75pc_epiboly[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="75pc_epiboly")])) ]         <- "Yes"
TE_GTF_sub_EXT$EXT_1_4_somites[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="1_4_somites")])) ]           <- "Yes"
TE_GTF_sub_EXT$EXT_14_19_somites[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="14_19_somites")])) ]       <- "Yes"
TE_GTF_sub_EXT$EXT_20_25_somites[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="20_25_somites")])) ]       <- "Yes"
TE_GTF_sub_EXT$EXT_prim_5[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="prim_5")])) ]                     <- "Yes"
TE_GTF_sub_EXT$EXT_prim_15[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="prim_15")])) ]                   <- "Yes"
TE_GTF_sub_EXT$EXT_prim_25[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="prim_25")])) ]                   <- "Yes"
TE_GTF_sub_EXT$EXT_long_pec[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="long_pec")])) ]                 <- "Yes"
TE_GTF_sub_EXT$EXT_protruding_mouth[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="protruding_mouth")])) ] <- "Yes"
TE_GTF_sub_EXT$EXT_day_4[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="day_4")])) ]                       <- "Yes"
TE_GTF_sub_EXT$EXT_day_5[ queryHits(findOverlaps(TE_GTF_sub_EXT, EXT_stages_noMT_CPM_0[which(names(EXT_stages_noMT_CPM_0)=="day_5")])) ]                       <- "Yes"

### Add Is_in_ext_3UTR column to TE_GTF_sub
TE_GTF_sub$Is_in_ext_3UTR <- rep("No",length(TE_GTF_sub_EXT))
TEs_ext_3UTR_idx <- unique(c(which(TE_GTF_sub_EXT$EXT_1_cell=='Yes'),which(TE_GTF_sub_EXT$EXT_2_cell=='Yes'),which(TE_GTF_sub_EXT$EXT_128_cell=='Yes'),which(TE_GTF_sub_EXT$EXT_1k_cell=='Yes'),which(TE_GTF_sub_EXT$EXT_dome=='Yes'),which(TE_GTF_sub_EXT$EXT_50pc_epiboly=='Yes'),which(TE_GTF_sub_EXT$EXT_shield=='Yes'),which(TE_GTF_sub_EXT$EXT_75pc_epiboly=='Yes'),which(TE_GTF_sub_EXT$EXT_1_4_somites=='Yes'),which(TE_GTF_sub_EXT$EXT_14_19_somites=='Yes'),which(TE_GTF_sub_EXT$EXT_20_25_somites=='Yes'),which(TE_GTF_sub_EXT$EXT_prim_5=='Yes'),which(TE_GTF_sub_EXT$EXT_prim_15=='Yes'),which(TE_GTF_sub_EXT$EXT_prim_25=='Yes'),which(TE_GTF_sub_EXT$EXT_long_pec=='Yes'),which(TE_GTF_sub_EXT$EXT_protruding_mouth=='Yes'),which(TE_GTF_sub_EXT$EXT_day_4=='Yes'),which(TE_GTF_sub_EXT$EXT_day_5=='Yes')))
TEs_ext_3UTR_idx <- TEs_ext_3UTR_idx[order(TEs_ext_3UTR_idx)]
TE_GTF_sub$Is_in_ext_3UTR[TEs_ext_3UTR_idx] <- "Yes"

### Save as GTF
export(TE_GTF_sub, opt$output)

#################