#!/bin/bash
# 4.run_DE_analysis.sh

### Run DESeq2 on Telescope counts together with gene counts
# Since gene expression accounts for a bigger fraction of the transcriptome, running DESeq 
# on TE counts together with gene counts ensures a better dispersion estimation, that will
# impact DESeq normalization
mkdir ./data/DE_analysis
Rscript ./scripts/run_DESeq2_Telescope.R

### DE analysis
Rscript ./scripts/run_DE_analysis_Telescope.R

###############