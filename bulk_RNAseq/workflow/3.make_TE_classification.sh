#!/bin/bash
# 3.make_TE_classification.sh

# Steps
#  - Get defrag file and parse it as necessary. 
#  - Run stringtie to get extended 3' UTRs
#    - Merge BAMs
#  - Run read-through categorization of TEs (R, ChipSeeker, ...)
#    - In the same script incorporate the extended 3' UTRs
#  - Merge telescope counts based on defrag TEs.

one_code_to_find_them_all.pl

### Get de-fragmented annotations
# To reconstructs fragmented repeats and full-length LTR elements from TE annotations we used the
# tool onecodetofindthemall.pl from Bailly-Bechet et al., 2014. DOI: https://doi.org/10.1186/1759-8753-5-13
# The output of this tool should be named danRer11.nonalt.fa.out.elem_sorted.csv
# and placed in the ./data direcotry.

# We then clean the output of onecodetofindthemall.pl to discard large TEs that have been joined
# based on the reference TE length.

time python clean_OCTFTA_output.py \
  --te_gtf /home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.gtf \
  --input /home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/defragmented/danRer11.nonalt.fa.out.elem_sorted.csv \
  --te_length /home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/defragmented/tes.lengths \
  --te_length_cutoff 2 \
  --output /home/qrovira/Projects/ZF_GRCz11_TEs/results/20_08_18_TEdefrag_overlap_DET/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.gtf
# Transform to BED
time python make_BED_from_clean_OCTFTA_GTF.py \
  --gtf /home/qrovira/Projects/ZF_GRCz11_TEs/results/20_08_18_TEdefrag_overlap_DET/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.gtf \
  --bed /home/qrovira/Projects/ZF_GRCz11_TEs/results/20_08_18_TEdefrag_overlap_DET/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.bed


### Run stringtie
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