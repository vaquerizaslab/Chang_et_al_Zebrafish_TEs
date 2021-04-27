#!/usr/bin/env Rscript
#run_DE_analysis_Telescope.R
# Perform DE analysis of TE loci with the count data from telescope.

suppressMessages(library(tictoc))
suppressMessages(library(DESeq2))

#######################
####  DE analysis  ####
# Read the results saved in: dds_TEloci_defrag.rds
tic("Read dds_TEloci_defrag.rds")
dds_TEloci_defrag <- readRDS("./data/DE_analysis/dds_TEloci_defrag.rds")
toc()

## Pairwise comparison (153)
# Send the pairwise comparison as jobs to run them in parallel
# The results of this pairwise comparison are NOT corrected for 
# multiple testing by the number of different comparisons.
coldata <- colData(dds_TEloci_defrag)
stages <- unique(coldata$Stage)
PAIRWISE_OUT_DIR <- "./data/DE_analysis/pairwise_results/non_multTestCorr"
dir.create("./data/DE_analysis/pairwise_results", showWarnings=FALSE)
dir.create(PAIRWISE_OUT_DIR, showWarnings=FALSE)
RES_SCRIPT <- "./scripts/deseq_results.R"
for (i in seq(1,length(stages))) {
        for (j in seq(1,length(stages))) {
                if(j > i){
                        writeLines(paste0("Rscript ",RES_SCRIPT," --dds ./data/DE_analysis/dds_TEloci_defrag.rds --cond Stage --cond_A ",stages[j]," --cond_B ",stages[i]," --output ",file.path(PAIRWISE_OUT_DIR,paste0("res.",stages[i],".",stages[j],".rds"))))
                        system(paste0("Rscript ",RES_SCRIPT," --dds ./data/DE_analysis/dds_TEloci_defrag.rds --cond Stage --cond_A ",stages[j]," --cond_B ",stages[i]," --output ",file.path(PAIRWISE_OUT_DIR,paste0("res.",stages[i],".",stages[j],".rds"))))
                }
        }
}

# Read the pairwise comparison output files
coldata <- colData(dds_TEloci_defrag)
stages <- unique(coldata$Stage)
#
pairwise_res_list <- list()
pairwise_names <- c()
counter <- 1
#
for (i in seq(1,length(stages))) {
        for (j in seq(1,length(stages))) {
                if(j > i){
                        writeLines(paste0("Reading: ",file.path(PAIRWISE_OUT_DIR,paste0("res.",stages[i],".",stages[j],".rds"))))
                        res_tmp <- readRDS( file.path(PAIRWISE_OUT_DIR,paste0("res.",stages[i],".",stages[j],".rds")) )
                        #
                        pairwise_res_list[[counter]] <- res_tmp
                        pairwise_names <- c(pairwise_names, paste0(stages[i],".",stages[j]))
                        counter <- counter+1
                }
        }
}
names(pairwise_res_list) <- pairwise_names

# Multiple testing correction
# Append all the results tables from all the comparisons one after the other, 
# and then use p.adjust() to calculate the adjusted p-value.
#
#Join result data frames into a single data frame
for (i in seq(1,length(pairwise_res_list))) {
        writeLines(paste0(i,"/153"))
        if (i == 1) {
                pairwise_res_DF <- pairwise_res_list[[i]]
                pairwise_res_DF$Comparison <- rep(pairwise_names[i],nrow(pairwise_res_DF))
        } else {
                pairwise_res_DF_tmp <- pairwise_res_list[[i]]
                pairwise_res_DF_tmp$Comparison <- rep(pairwise_names[i],nrow(pairwise_res_DF_tmp))
                pairwise_res_DF <- rbind(pairwise_res_DF, pairwise_res_DF_tmp)
        }
}
#Re-calculate p-adjusted value
pairwise_res_DF$padj_2 <- p.adjust(pairwise_res_DF$pvalue, method="BH")
#Save re-calculated p-adjusted value in every data frame
for (i in seq(1,length(pairwise_res_list))) {
        writeLines(paste0(i,"/153"))
        pairwise_res_list[[i]]$padj_all <- pairwise_res_DF$padj_2[which(pairwise_res_DF$Comparison==pairwise_names[i])]
}
#Save result table as .rds again
PAIRWISE_OUT_DIR <- "./data/DE_analysis/pairwise_results/multTestCorr"
dir.create(PAIRWISE_OUT_DIR, showWarnings=FALSE)
for (i in seq(1,length(stages))) {
        for (j in seq(1,length(stages))) {
                if (paste0(stages[i],".",stages[j]) %in% pairwise_names) {
                        idx <- which(pairwise_names==paste0(stages[i],".",stages[j]))
                        writeLines(paste0(idx,"/153"))
                        res_tmp <- pairwise_res_list[[idx]]
                        saveRDS(res_tmp, file.path(PAIRWISE_OUT_DIR,paste0("res.",stages[i],".",stages[j],".rds")))
                }
        }
}


### Get the significantly DE TEs in any pairwise comparison
stages <- c("1_cell","2_cell","128_cell","1k_cell","dome","50pc_epiboly","shield","75pc_epiboly","1_4_somites","14_19_somites","20_25_somites","prim_5","prim_15","prim_25","long_pec","protruding_mouth","day_4","day_5")
DE_TE_DIR <- "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/DE_analysis_TEdefrag/DE_TE_loci"
SIG_DE_DIR <- "./data/DE_analysis/pairwise_results/sig_DE"
dir.create(SIG_DE_DIR, showWarnings=FALSE)
#
padj <- 0.01
log2FC <- 0
writeLines(paste0("p-adjusted<",padj," & log2FC>=",log2FC))
tic()
counter <- 1
list_sigDE_TEloci <- list()
for (i in seq(1,length(stages))) {
        for (j in seq(1,length(stages))) {
                if(j > i){
                        #writeLines(paste0("Reading: res.",stages[i],".",stages[j]))
                        res_Teles_DE_TEloci   <- readRDS(file.path("./data/DE_analysis/pairwise_results/multTestCorr", paste0("res.",stages[i],".",stages[j],".rds")))
                        sigDE_Teles_DE_TEloci <- rownames(res_Teles_DE_TEloci[which(res_Teles_DE_TEloci$padj_all<padj & abs(res_Teles_DE_TEloci$log2FoldChange)>=log2FC),])
                        list_sigDE_TEloci[[counter]] <- sigDE_Teles_DE_TEloci
                        counter <- counter+1
                }
        }
}
#
sigDE_TEloci_anyPair <- unique(unlist(list_sigDE_TEloci))
write(sigDE_TEloci_anyPair, file.path("./data/DE_analysis/pairwise_results/sig_DE", paste0("All_sigDE_TEloci.pAdj",padj,".log2FC",log2FC,".txt")))


#############