#!/bin/bash
# 3.make_TE_classification.sh

# Steps
#  - Get defrag file and parse it as necessary. 
#  - Run stringtie to get extended 3' UTRs
#    - Merge BAMs
#  - Run read-through categorization of TEs (R, ChipSeeker, ...)
#    - In the same script incorporate the extended 3' UTRs
#  - Merge telescope counts based on defrag TEs.


### Get de-fragmented annotations
# To reconstructs fragmented repeats and full-length LTR elements from TE annotations we used the
# tool onecodetofindthemall.pl from Bailly-Bechet et al., 2014. DOI: https://doi.org/10.1186/1759-8753-5-13
# The output of this tool should be named danRer11.nonalt.fa.out.elem_sorted.csv
# and placed in the ./data direcotry.
# A 2 column file with the length of the consensus sequence for each TE is also necessary: ./data/tes.lengths

# We then clean the output of onecodetofindthemall.pl to discard large TEs that have been joined
# based on the reference TE length.

TE_GTF="./data/danRer11.TEtrans_uID.gtf"
time python ./scripts/clean_OCTFTA_output.py \
  --te_gtf ./data/danRer11.TEtrans_uID.gtf \
  --input ./data/danRer11.nonalt.fa.out.elem_sorted.csv \
  --te_length ./data/tes.lengths \
  --te_length_cutoff 2 \
  --output ./data/danRer11.TEtrans_uID.dfrag.gtf
# Transform to BED
time python ./scripts/make_BED_from_clean_OCTFTA_GTF.py \
  --gtf ./data/danRer11.TEtrans_uID.dfrag.gtf \
  --bed ./data/danRer11.TEtrans_uID.dfrag.bed



### Run stringtie
# In order to have more sequencing depth, replicates were merged before transcriptome assembly.
# Merge BAMs


# Run Stringtie (only with the good parameters)
# /home/qrovira/Projects/ZF_GRCz11_TEs/results/20_05_18_White2017_Stringtie/run_Stringtie.sh
# Get extended 3' UTRs
# get_pervasive_tails.R


### Make TE categories 
# ChipSeeker
# Add get_pervasive_tails results
# Export as new GTF

### Merge telescope counts
#/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/join_TEcounts_defragID.sh
#/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/merge_Teles_counts/merge_Teles_counts_defrag.py

##################