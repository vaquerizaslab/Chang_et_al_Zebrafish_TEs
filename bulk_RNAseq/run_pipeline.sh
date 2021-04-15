#!/bin/bash
# run_pipeline.sh

# This script contains the step by step parts of the bulk RNA-seq analysis.
# This script is not meant to be executed and reproduce all the resuts, 
# instead, each step is supposed to be run independently.

### 1. Fastq trimming and mapping
bash workflow/1.fastq_trim_and_map.sh

### 2. Gene and TE counts
bash workflow/2.make_gene_and_TE_counts.sh

### 3. TE classification
bash workflow/3.make_TE_classification.sh
#(include the merge defrag step as well)

### 4. DE analysis
bash workflow/4.run_DE_analysis.sh

### 5. Clustering
bash workflow/5.TE_clustering.sh

### 6. Enrichment
bash workflow/6.TE_enrichment.sh

### 7. Generate figure plots
bash workflow/7.figure_plots.sh

####################