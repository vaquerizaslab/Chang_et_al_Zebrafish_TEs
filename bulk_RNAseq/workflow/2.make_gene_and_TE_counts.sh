#!/bin/bash
# 2.make_gene_and_TE_counts.sh

# Steps
#  - Go from supp-data-1/danRer11.nonalt.fa.out to the TEtranscript ready GTF
#  - Run TE transcripts
#  - Run Telescope
#  - XXX


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


### Get gene annotations
# Download GRCz11.98 gene annotations from Ensembl
wget -O ./data/Danio_rerio.GRCz11.98.chr.gtf.gz http://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.chr.gtf.gz
# Download ERCC
cd data
wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip
rm ERCC92.zip ERCC92.fa
cd ..
# Add ERCC gtf to GRCz11.98 annotations
cat <(zcat ./data/Danio_rerio.GRCz11.98.chr.gtf.gz) \
    ./data/ERCC92.gtf \
    | gzip > ./data/Danio_rerio.GRCz11.98.chr.ERCC.gtf.gz
# BAM files chromosomes have the "chr" on it while the annotations don't.
# Add "chr" to each entry of the GTF but the ERCCs. Take special care of MT chromosome.
zcat ./data/Danio_rerio.GRCz11.98.chr.ERCC.gtf.gz \
  | sed -e 's/^/chr/' \
  | sed -e 's/^chr#/#/;s/^chrERCC/ERCC/;s/^chrMT/chrM/' \
  | gzip > ./data/Danio_rerio.GRCz11.98.chr.ERCC.addChr.gtf.gz


### Run TEtranscripts
# Run samples in pairs
BAM_DIR="./data/BAMs"
GE_GTF="./data/Danio_rerio.GRCz11.98.chr.ERCC.addChr.gtf.gz"
TE_GTF="./data/danRer11.TEtrans_uID.gtf"
TETRANS_OUT_DIR="./data/TEtranscripts"
mkdir $TETRANS_OUT_DIR
# 1_cell        vs      2_cell
S1="1_cell"; S2="2_cell"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# 128_cell      vs      1k_cell
S1="128_cell"; S2="1k_cell"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# dome          vs      50pc_epiboly
S1="dome"; S2="50pc_epiboly"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# shield        vs      75pc_epiboly
S1="shield"; S2="75pc_epiboly"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# 1_4_somites   vs      14_19_somites
S1="1_4_somites"; S2="14_19_somites"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# 20_25_somites vs      prim_5
S1="20_25_somites"; S2="prim_5"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# prim_15       vs      prim_25
S1="prim_15"; S2="prim_25"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# long_pec      vs      protruding_mouth
S1="long_pec"; S2="protruding_mouth"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# day_4         vs      day_5
S1="day_4"; S2="day_5"
TEtranscripts --format BAM --stranded reverse --mode multi --sortByPos --control $BAM_DIR/${S1}_rep_1.bam $BAM_DIR/${S1}_rep_2.bam $BAM_DIR/${S1}_rep_3.bam $BAM_DIR/${S1}_rep_4.bam $BAM_DIR/${S1}_rep_5.bam --treatment $BAM_DIR/${S2}_rep_1.bam $BAM_DIR/${S2}_rep_2.bam $BAM_DIR/${S2}_rep_3.bam $BAM_DIR/${S2}_rep_4.bam $BAM_DIR/${S2}_rep_5.bam --GTF $GE_GTF --TE $TE_GTF --project ${TETRANS_OUT_DIR}/${S1}_vs_${S2}
# Join the *.cntTable files
DIR_PATH="\.\/data\/BAMs\/"
paste <(cut -f1 ./data/TEtranscripts/1_cell_vs_2_cell.cntTable | sed 's/gene\/TE/Feature/') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/1_cell_vs_2_cell.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/128_cell_vs_1k_cell.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/dome_vs_50pc_epiboly.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/shield_vs_75pc_epiboly.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/1_4_somites_vs_14_19_somites.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/20_25_somites_vs_prim_5.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/prim_15_vs_prim_25.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/long_pec_vs_protruding_mouth.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
      <(sed "s/${DIR_PATH}//g" ./data/TEtranscripts/day_4_vs_day_5.cntTable | sed -e 's/.bam.[TC]//g' | sed 's/gene\/TE/Feature/'| awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}') \
  | sed 's/"//g' \
  > ./data/TEtranscripts/join.cntTable

# Load join.cntTable into R as DESeq object and save it in RDS format
R ./scripts/load_TEtrans_cntTable.R


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