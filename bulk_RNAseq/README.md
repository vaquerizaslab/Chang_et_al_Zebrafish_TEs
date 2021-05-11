# Bulk RNA-seq analysis code for Chang et al. 2021

This analysis includes the following steps:
* Mapp reads
* Count gene and TEs
* TE classification
* TE DE analysis
* Clustering
* TE Enrichment

### Structure
The code should be reproduced following the steps on the `run_pipeline.sh` script, which calls the scripts for each specific step of the analysis. The rest of the scripts are save in the `scripts/` directory. 
The pipeline assumes that the `data/` directory contains a folder `fastq/` with the fastq files from White et al. 2017 with the following sample names:
```
1_cell
2_cell
128_cell
1k_cell
dome
50pc_epiboly
shield
75pc_epiboly
1_4_somites
14_19_somites
20_25_somites
prim_5
prim_15
prim_25
long_pec
protruding_mouth
day_4
day_5
```
and replicate suffix:
```
_rep_1_R1.fastq.gz
_rep_1_R2.fastq.gz
_rep_2_R1.fastq.gz
_rep_2_R2.fastq.gz
_rep_3_R1.fastq.gz
_rep_3_R2.fastq.gz
_rep_4_R1.fastq.gz
_rep_4_R2.fastq.gz
_rep_5_R1.fastq.gz
_rep_5_R2.fastq.gz
```
i.e. `dome_rep_3_R1.fastq.gz` & `dome_rep_3_R2.fastq.gz`

Large publicly available files are not included and should be downloaded independently.