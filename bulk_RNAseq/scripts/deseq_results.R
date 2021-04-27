#!/usr/bin/env Rscript
# deseq_results.R

suppressMessages(library("optparse"))
option_list <- list(
  make_option(c("--dds"), type="character", default=NULL, 
              help="Path to DESeqDataSet file in .rds format"),
  make_option(c("--cond"), type="character", default=NULL, 
              help="Name of the condition column to be used (it has to match the colData column name). It will run like:\n
              results(dds, contrast=c(--cond ,--cond_A ,--cond_B))"),
  make_option(c("--cond_A"), type="character", default=NULL, 
              help="Name of the condition A to be used. It will run like:\n
              results(dds, contrast=c(--cond ,--cond_A ,--cond_B))"),
  make_option(c("--cond_B"), type="character", default=NULL, 
              help="Name of the condition A to be used. It will run like:\n
              results(dds, contrast=c(--cond ,--cond_A ,--cond_B))"),
  make_option(c("--output"), type="character", default=NULL, 
              help="Path to output file in .rds")
);
opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

#### Print parameters
writeLines( paste0("Parameter --dds: ", opt$dds) )
writeLines( paste0("Parameter --cond: ", opt$cond) )
writeLines( paste0("Parameter --cond_A: ", opt$cond_A) )
writeLines( paste0("Parameter --cond_B: ", opt$cond_B) )
writeLines( paste0("Parameter --output: ", opt$output) )

#### Loading packages
suppressMessages(library(tictoc))
tic("Loading packages")
suppressMessages(library(DESeq2))
toc()

#### Reading dds file
tic(paste0("Reading file: ",opt$dds))
dds <- readRDS(opt$dds)
toc()

#### Running DESeq2 results()
tic(paste0("Results ",opt$cond," ",opt$cond_A," vs ",opt$cond_B))
res <- results(dds, contrast=c(opt$cond, opt$cond_A, opt$cond_B))
toc()

#### Saving results()
tic(paste0("Saving results: ",opt$output))
saveRDS(res, opt$output)
toc()

############