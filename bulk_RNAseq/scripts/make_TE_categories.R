#!/usr/bin/env Rscript
#make_TE_categories.R
# Make TE classification
# Do NOT consider if a TE is DE or not for it.
# Make 2 classifications:
#  - Detailed:
#    + Ex:        Exon (SS)
#    + In_g:      Intron (SS) gene expression
#    + In_ng:     Intron (SS) no gene expression
#    + 3utr_ext:  Pervasive transcription (non-annotated 3' UTRs)
#    + nRT:       Non read-through
#  - Non detailed:
#    + Gene_dependent: Ex + In_g + 3utr_ext
#    + Self_expressed: In_ng + nRT
#
# What to do when different parts of an assembled
# TE are in different categories?
#   Give all the parts the most restrictive classification.
#     Ex > 3utr_ext > In_g > In_ng > nRT
#     Gene_dependent > Self_expressed


suppressWarnings(suppressMessages(library(optparse)))
option_list <- list(
  make_option(c("-t","--te_gtf"), type="character", default=NULL, 
              help="Path to TE GTF."),
  make_option(c("-g","--te_gtf_ext"), type="character", default=NULL, 
              help="Path to TE GTF with extended 3' UTR info."),
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

### Load TE annotations with extended 3' UTR data
tic("Loading TE annotations with extended 3' UTR data")
TE_GTF_defrag_ext <- import(opt$te_gtf_ext)
toc(); writeLines("")

### Add Is_in_ext_3UTR column to TE_GTF_defrag
TE_GTF_defrag$Is_in_ext_3UTR <- TE_GTF_defrag_ext$Is_in_ext_3UTR[match(TE_GTF_defrag$transcript_id, TE_GTF_defrag_ext$transcript_id)]
rm(TE_GTF_defrag_ext)

### Load Gene expression matrix
tic("Read ./data/TEtranscripts/dds_TEtrans_GE.rds")
dds_TEtrans_GE <- readRDS("./data/TEtranscripts/dds_TEtrans_GE.rds")
gene_counts <- counts(dds_TEtrans_GE, normalized=TRUE)
#
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
gene_counts_meanRep <- mean_by_rep(gene_counts)
toc(); writeLines("")

### Subset expressed genes
# Consider a gene expressed if it has more than 10 counts in at least 1 stage
tic("Subset expressed genes")
expressed_genes <- rownames(gene_counts_meanRep)[!(apply(gene_counts_meanRep,1,function(x) all(x<10)))]
toc(); writeLines("")

### Make read-through categories
tic("Make read-through categories")
#   Ex:         Exon (SS)
#   In_g:       Intron (SS) gene expression
#   In_ng:      Intron (SS) no gene expression
#   3utr_ext:   Non-annotated 3' UTRs
#   nRT:        Non read-through
#   Ex > 3utr_ext > In_g > In_ng > nRT
# Create following columns:
#  - Is exon overlapping
#  - Is intron overlapping
#  - Is intron gene expressed
#  - Is in extended 3' UTR (already exist)
# If NO in all, then it's considered "Non read-through"
TE_GTF_defrag_tmp <- TE_GTF_defrag
TE_GTF_defrag_tmp$Is_EX_overl <- "No"
TE_GTF_defrag_tmp$Is_EX_overl[which(TE_GTF_defrag_tmp$annotation_3=="Exon.SS" |
                                    TE_GTF_defrag_tmp$annotation_3=="3_UTR.SS" |
                                    TE_GTF_defrag_tmp$annotation_3=="5_UTR.SS")] <- "Yes"
#
TE_GTF_defrag_tmp$Is_IN_overl <- "No"
TE_GTF_defrag_tmp$Is_IN_overl[which(TE_GTF_defrag_tmp$annotation_3=="Intron.SS")] <- "Yes"
#
TE_GTF_defrag_tmp$Is_IN_overl_gene_expr <- "No"
TE_GTF_defrag_tmp$Is_IN_overl_gene_expr[which(TE_GTF_defrag_tmp$Is_IN_overl=="Yes" & TE_GTF_defrag_tmp$geneId %in% expressed_genes)] <- "Yes"
#table(TE_GTF_defrag_tmp$Is_IN_overl,TE_GTF_defrag_tmp$Is_IN_overl_gene_expr)
#           No     Yes
#  No  2794615       0
#  Yes   97497  908323
#
# Ex:        Is_EX_overl=="Yes"
# 3utr_ext:  Is_EX_overl=="No"  & Is_in_ext_3UTR=="Yes"
# In_g:      Is_EX_overl=="No"  & Is_in_ext_3UTR=="No"  & Is_IN_overl_gene_expr=="Yes"
# In_ng:     Is_EX_overl=="No"  & Is_in_ext_3UTR=="No"  & Is_IN_overl_gene_expr=="No"  & Is_IN_overl=="Yes"
# nRT:       Is_EX_overl=="No"  & Is_in_ext_3UTR=="No"  & Is_IN_overl_gene_expr=="No"  & Is_IN_overl=="No"
TE_GTF_defrag_tmp$Rt_catg_detail <- "BLANK"
TE_GTF_defrag_tmp$Rt_catg_detail[which(TE_GTF_defrag_tmp$Is_EX_overl=="Yes")] <- "Ex"
TE_GTF_defrag_tmp$Rt_catg_detail[which(TE_GTF_defrag_tmp$Is_EX_overl=="No"  & TE_GTF_defrag_tmp$Is_in_ext_3UTR=="Yes")] <- "3utr_ext"
TE_GTF_defrag_tmp$Rt_catg_detail[which(TE_GTF_defrag_tmp$Is_EX_overl=="No"  & TE_GTF_defrag_tmp$Is_in_ext_3UTR=="No"  & TE_GTF_defrag_tmp$Is_IN_overl_gene_expr=="Yes")] <- "In_g"
TE_GTF_defrag_tmp$Rt_catg_detail[which(TE_GTF_defrag_tmp$Is_EX_overl=="No"  & TE_GTF_defrag_tmp$Is_in_ext_3UTR=="No"  & TE_GTF_defrag_tmp$Is_IN_overl_gene_expr=="No"  & TE_GTF_defrag_tmp$Is_IN_overl=="Yes")] <- "In_ng"
TE_GTF_defrag_tmp$Rt_catg_detail[which(TE_GTF_defrag_tmp$Is_EX_overl=="No"  & TE_GTF_defrag_tmp$Is_in_ext_3UTR=="No"  & TE_GTF_defrag_tmp$Is_IN_overl_gene_expr=="No"  & TE_GTF_defrag_tmp$Is_IN_overl=="No")] <- "nRT"
#
TE_GTF_defrag_tmp$Rt_catg <- "BLANK"
TE_GTF_defrag_tmp$Rt_catg[which(TE_GTF_defrag_tmp$Rt_catg_detail=="Ex" | 
                                TE_GTF_defrag_tmp$Rt_catg_detail=="3utr_ext" | 
                                TE_GTF_defrag_tmp$Rt_catg_detail=="In_g")] <- "Gene_dependent"
TE_GTF_defrag_tmp$Rt_catg[which(TE_GTF_defrag_tmp$Rt_catg_detail=="In_ng" | 
                                TE_GTF_defrag_tmp$Rt_catg_detail=="nRT")] <- "Self_expressed"
#
TE_GTF_defrag <- TE_GTF_defrag[,1:15]
TE_GTF_defrag$Rt_catg        <- TE_GTF_defrag_tmp$Rt_catg
TE_GTF_defrag$Rt_catg_detail <- TE_GTF_defrag_tmp$Rt_catg_detail
toc(); writeLines("")

### Fix read-through categories for defrag cases
TE_GTF_defrag_DF <- as.data.frame(TE_GTF_defrag)
tic("Dplyr concatenate RT by defrag_id")
TE_GTF_defrag_DF_c <- TE_GTF_defrag_DF %>%
    group_by(defrag_transcript_id) %>%
    dplyr::mutate(Rt_catg_c=paste0(Rt_catg, collapse=":"),
                  Rt_catg_detail_c=paste0(Rt_catg_detail, collapse=":"))
toc()
#
Gene_or_Self_expressed <- function(vector){
	if("Gene_dependent" %in% vector){
		out_string <- "Gene_dependent"
	} else{
		out_string <- "Self_expressed"
	}
	return(out_string)
}
Detailed_RT_hierarchy <- function(vector){
	if(length(vector)==1){
		out_string <- vector[1]
	} else{
		if("Ex" %in% vector){
			out_string <- "Ex"
		} else if("3utr_ext" %in% vector){
			out_string <- "3utr_ext"
		} else if("In_g" %in% vector){
			out_string <- "In_g"
		} else if("In_ng" %in% vector){
			out_string <- "In_ng"
		} else{
			out_string <- "nRT"
		}
	}
	return(out_string)
}
#
TE_GTF_defrag_DF_c$Rt_catg_c2 <- unlist(lapply(strsplit(TE_GTF_defrag_DF_c$Rt_catg_c,":"), Gene_or_Self_expressed))
TE_GTF_defrag_DF_c$Rt_catg_detail_c2 <- unlist(lapply(strsplit(TE_GTF_defrag_DF_c$Rt_catg_detail_c,":"), Detailed_RT_hierarchy))
#
TE_GTF_defrag$Rt_catg_noDefrag <- TE_GTF_defrag$Rt_catg
TE_GTF_defrag$Rt_catg_detail_noDefrag <- TE_GTF_defrag$Rt_catg_detail
TE_GTF_defrag$Rt_catg        <- TE_GTF_defrag_DF_c$Rt_catg_c2
TE_GTF_defrag$Rt_catg_detail <- TE_GTF_defrag_DF_c$Rt_catg_detail_c2

### Save as GTF and RDS
tic("Save as GTF and RDS")
saveRDS(TE_GTF_defrag, sub(".gtf",".rds",opt$output))
export(TE_GTF_defrag, opt$output)
toc(); writeLines("")

##############