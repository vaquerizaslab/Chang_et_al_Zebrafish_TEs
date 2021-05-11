# Analysis code for Chang et al. 2021

This repository contains the necessary code to reproduce the analysis performed in Chang et al. 2021, '[A genomic portrait of zebrafish transposable elements and their spatiotemporal embryonic expression](https://www.biorxiv.org/content/10.1101/2021.04.08.439009v1.full)'.

The code provided here complements and expands the Methods described in the paper with the aim to increase the reproducibility of our work and help understand the analysis performed. Note that this repository contains the code, and is not a container to execute the entire workflow af the paper.

## Overview
The analysis are split between three separate sections:
 * Genomic TE analysis
 * Bulk RNA-seq analysis
 * scRNA-seq analysis.  

These sections were developed independently, but relying on the same genome version and gene and TE annotations.


### Data
In this work we analyzed public available bulk RNA-seq data from [White et al. 2017](https://elifesciences.org/articles/30860) with ENA accession ERP014517, and scRNA-seq data from [Farrell et al. 2019](https://science.sciencemag.org/content/360/6392/eaar3131) with GEO accession GSE106587.


### Software requirements
Different sections might have different software requirements. Please check each section.


### Genomic TE analysis
The code to reproduce the genomic TE analysis can be found in the directory `genomic_TE`.


### Bulk RNA-seq
The code to reproduce the bulk RNA-seq analysis can be found in the directory `bulk_RNAseq`.


### scRNA-seq
The code to reproduce the scRNA-seq analysis can be found in the directory `scRNAseq`.


### Contact
If you have any question regarding the code, please open an issue or contact us.
