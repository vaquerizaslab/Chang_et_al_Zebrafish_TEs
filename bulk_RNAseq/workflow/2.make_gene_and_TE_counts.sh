#!/bin/bash
# 2.make_gene_and_TE_counts.sh

# Steps
#  - Go from supp-data-1/danRer11.nonalt.fa.out to the TEtranscript ready GTF
#  - Run TE transcripts
#  - Run Telescope
#  - XXX


#lh /home/qrovira/Projects/ZF_GRCz11_TEs/results/20_01_07_White2017_TEtranscripts_geneAnnot
#-rw-r--r-- 1 qrovira qrovira 3.4K Jan  7  2020 adapt_GTF_TEtranscript_fix_TEfam.py
#-rw-r--r-- 1 qrovira qrovira  18K Jan  7  2020 run_White2017_TEtranscripts_geneAnnot.sh
#copy danRer11_classes_wredundant.txt to data/


### Prepare TE annotations
# Get the RepeatMaker output file from Supplementary Data 1 and store it ./data
# Use rmsk2bed script from BEDOPS (version 2.4.27) to convert .fa.out to BED file
rmsk2bed < ./data/danRer11.nonalt.fa.out > ./data/danRer11.nonalt.fa.out.bed
# Convert to GTF format
time python ./scripts/RM_bed2GTF.py -i ./data/danRer11.nonalt.fa.out.bed -o ./data/danRer11.nonalt.fa.out.gtf
# Adapt GTF to TEtranscript requested format
# Use ./data/danRer11_classes_wredundant.txt file to solve some TE superfamily/family name redundancy.
# Subset TE classes: DNA, LINE, LTR, RC, SINE
sed 's/?//g' ./data/danRer11_classes_wredundant.txt \
  | awk '{if ($2=="DNA" || $2=="LINE" || $2=="LTR" || $2=="RC" || $2=="SINE") print }' \
  > ./data/danRer11_classes_wredundant.subClass.txt
python ./scripts/adapt_GTF_TEtranscript.py \
  -i ./data/danRer11.nonalt.fa.out.gtf \
  -c ./data/danRer11_classes_wredundant.subClass.txt \
  -o ./data/danRer11.TEtrans_uID.gtf


### Run TEtranscripts
# Run pairwise

# Join the *.cntTable files

# Load them into R and save a DESeq object as RDS?


### Run Telescope
# Split BAMs
# 20_04_23_White2017_Split_BAM_by_strand/run_White2017_Split_BAM_by_strand.sh

# Run telescope on split BAMs
# 20_04_23_White2017_Telescope_stranded/run_White2017_Telescope_stranded.sh

# Sort all files on the same way

# Join results into a single table


# Merge fwd/rev results? (I think I did this step in the DE analysis part)
# Load them into R and save a DESeq object as RDS?




###########