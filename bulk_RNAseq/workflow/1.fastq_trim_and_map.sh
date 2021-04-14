#!/bin/bash
# 1.fastq_trim_and_map.sh

### Fastq files from White et al. 2017 should be in data/fastq/ with the appropriate name.
ls data/fastq/*
#128_cell_rep_1_R1.fastq.gz       1k_cell_rep_3_R2.fastq.gz        day_4_rep_1_R1.fastq.gz     prim_15_rep_3_R2.fastq.gz
#128_cell_rep_1_R2.fastq.gz       1k_cell_rep_4_R1.fastq.gz        day_4_rep_1_R2.fastq.gz     prim_15_rep_4_R1.fastq.gz
#128_cell_rep_2_R1.fastq.gz       1k_cell_rep_4_R2.fastq.gz        day_4_rep_2_R1.fastq.gz     prim_15_rep_4_R2.fastq.gz
#128_cell_rep_2_R2.fastq.gz       1k_cell_rep_5_R1.fastq.gz        day_4_rep_2_R2.fastq.gz     prim_15_rep_5_R1.fastq.gz
#128_cell_rep_3_R1.fastq.gz       1k_cell_rep_5_R2.fastq.gz        day_4_rep_3_R1.fastq.gz     prim_15_rep_5_R2.fastq.gz
#128_cell_rep_3_R2.fastq.gz       20_25_somites_rep_1_R1.fastq.gz  day_4_rep_3_R2.fastq.gz     prim_25_rep_1_R1.fastq.gz
#128_cell_rep_4_R1.fastq.gz       20_25_somites_rep_1_R2.fastq.gz  day_4_rep_4_R1.fastq.gz     prim_25_rep_1_R2.fastq.gz
#128_cell_rep_4_R2.fastq.gz       20_25_somites_rep_2_R1.fastq.gz  day_4_rep_4_R2.fastq.gz     prim_25_rep_2_R1.fastq.gz
#128_cell_rep_5_R1.fastq.gz       20_25_somites_rep_2_R2.fastq.gz  day_4_rep_5_R1.fastq.gz     prim_25_rep_2_R2.fastq.gz
#128_cell_rep_5_R2.fastq.gz       20_25_somites_rep_3_R1.fastq.gz  day_4_rep_5_R2.fastq.gz     prim_25_rep_3_R1.fastq.gz
#14_19_somites_rep_1_R1.fastq.gz  20_25_somites_rep_3_R2.fastq.gz  day_5_rep_1_R1.fastq.gz     prim_25_rep_3_R2.fastq.gz
#14_19_somites_rep_1_R2.fastq.gz  20_25_somites_rep_4_R1.fastq.gz  day_5_rep_1_R2.fastq.gz     prim_25_rep_4_R1.fastq.gz
#14_19_somites_rep_2_R1.fastq.gz  20_25_somites_rep_4_R2.fastq.gz  day_5_rep_2_R1.fastq.gz     prim_25_rep_4_R2.fastq.gz
#14_19_somites_rep_2_R2.fastq.gz  20_25_somites_rep_5_R1.fastq.gz  day_5_rep_2_R2.fastq.gz     prim_25_rep_5_R1.fastq.gz
#14_19_somites_rep_3_R1.fastq.gz  20_25_somites_rep_5_R2.fastq.gz  day_5_rep_3_R1.fastq.gz     prim_25_rep_5_R2.fastq.gz
#14_19_somites_rep_3_R2.fastq.gz  2_cell_rep_1_R1.fastq.gz         day_5_rep_3_R2.fastq.gz     prim_5_rep_1_R1.fastq.gz
#14_19_somites_rep_4_R1.fastq.gz  2_cell_rep_1_R2.fastq.gz         day_5_rep_4_R1.fastq.gz     prim_5_rep_1_R2.fastq.gz
#14_19_somites_rep_4_R2.fastq.gz  2_cell_rep_2_R1.fastq.gz         day_5_rep_4_R2.fastq.gz     prim_5_rep_2_R1.fastq.gz
#14_19_somites_rep_5_R1.fastq.gz  2_cell_rep_2_R2.fastq.gz         day_5_rep_5_R1.fastq.gz     prim_5_rep_2_R2.fastq.gz
#14_19_somites_rep_5_R2.fastq.gz  2_cell_rep_3_R1.fastq.gz         day_5_rep_5_R2.fastq.gz     prim_5_rep_3_R1.fastq.gz
#1_4_somites_rep_1_R1.fastq.gz    2_cell_rep_3_R2.fastq.gz         dome_rep_1_R1.fastq.gz      prim_5_rep_3_R2.fastq.gz
#1_4_somites_rep_1_R2.fastq.gz    2_cell_rep_4_R1.fastq.gz         dome_rep_1_R2.fastq.gz      prim_5_rep_4_R1.fastq.gz
#1_4_somites_rep_2_R1.fastq.gz    2_cell_rep_4_R2.fastq.gz         dome_rep_2_R1.fastq.gz      prim_5_rep_4_R2.fastq.gz
#1_4_somites_rep_2_R2.fastq.gz    2_cell_rep_5_R1.fastq.gz         dome_rep_2_R2.fastq.gz      prim_5_rep_5_R1.fastq.gz
#1_4_somites_rep_3_R1.fastq.gz    2_cell_rep_5_R2.fastq.gz         dome_rep_3_R1.fastq.gz      prim_5_rep_5_R2.fastq.gz
#1_4_somites_rep_3_R2.fastq.gz    50pc_epiboly_rep_1_R1.fastq.gz   dome_rep_3_R2.fastq.gz      protruding_mouth_rep_1_R1.fastq.gz
#1_4_somites_rep_4_R1.fastq.gz    50pc_epiboly_rep_1_R2.fastq.gz   dome_rep_4_R1.fastq.gz      protruding_mouth_rep_1_R2.fastq.gz
#1_4_somites_rep_4_R2.fastq.gz    50pc_epiboly_rep_2_R1.fastq.gz   dome_rep_4_R2.fastq.gz      protruding_mouth_rep_2_R1.fastq.gz
#1_4_somites_rep_5_R1.fastq.gz    50pc_epiboly_rep_2_R2.fastq.gz   dome_rep_5_R1.fastq.gz      protruding_mouth_rep_2_R2.fastq.gz
#1_4_somites_rep_5_R2.fastq.gz    50pc_epiboly_rep_3_R1.fastq.gz   dome_rep_5_R2.fastq.gz      protruding_mouth_rep_3_R1.fastq.gz
#1_cell_rep_1_R1.fastq.gz         50pc_epiboly_rep_3_R2.fastq.gz   long_pec_rep_1_R1.fastq.gz  protruding_mouth_rep_3_R2.fastq.gz
#1_cell_rep_1_R2.fastq.gz         50pc_epiboly_rep_4_R1.fastq.gz   long_pec_rep_1_R2.fastq.gz  protruding_mouth_rep_4_R1.fastq.gz
#1_cell_rep_2_R1.fastq.gz         50pc_epiboly_rep_4_R2.fastq.gz   long_pec_rep_2_R1.fastq.gz  protruding_mouth_rep_4_R2.fastq.gz
#1_cell_rep_2_R2.fastq.gz         50pc_epiboly_rep_5_R1.fastq.gz   long_pec_rep_2_R2.fastq.gz  protruding_mouth_rep_5_R1.fastq.gz
#1_cell_rep_3_R1.fastq.gz         50pc_epiboly_rep_5_R2.fastq.gz   long_pec_rep_3_R1.fastq.gz  protruding_mouth_rep_5_R2.fastq.gz
#1_cell_rep_3_R2.fastq.gz         75pc_epiboly_rep_1_R1.fastq.gz   long_pec_rep_3_R2.fastq.gz  shield_rep_1_R1.fastq.gz
#1_cell_rep_4_R1.fastq.gz         75pc_epiboly_rep_1_R2.fastq.gz   long_pec_rep_4_R1.fastq.gz  shield_rep_1_R2.fastq.gz
#1_cell_rep_4_R2.fastq.gz         75pc_epiboly_rep_2_R1.fastq.gz   long_pec_rep_4_R2.fastq.gz  shield_rep_2_R1.fastq.gz
#1_cell_rep_5_R1.fastq.gz         75pc_epiboly_rep_2_R2.fastq.gz   long_pec_rep_5_R1.fastq.gz  shield_rep_2_R2.fastq.gz
#1_cell_rep_5_R2.fastq.gz         75pc_epiboly_rep_3_R1.fastq.gz   long_pec_rep_5_R2.fastq.gz  shield_rep_3_R1.fastq.gz
#1k_cell_rep_1_R1.fastq.gz        75pc_epiboly_rep_3_R2.fastq.gz   prim_15_rep_1_R1.fastq.gz   shield_rep_3_R2.fastq.gz
#1k_cell_rep_1_R2.fastq.gz        75pc_epiboly_rep_4_R1.fastq.gz   prim_15_rep_1_R2.fastq.gz   shield_rep_4_R1.fastq.gz
#1k_cell_rep_2_R1.fastq.gz        75pc_epiboly_rep_4_R2.fastq.gz   prim_15_rep_2_R1.fastq.gz   shield_rep_4_R2.fastq.gz
#1k_cell_rep_2_R2.fastq.gz        75pc_epiboly_rep_5_R1.fastq.gz   prim_15_rep_2_R2.fastq.gz   shield_rep_5_R1.fastq.gz
#1k_cell_rep_3_R1.fastq.gz        75pc_epiboly_rep_5_R2.fastq.gz   prim_15_rep_3_R1.fastq.gz   shield_rep_5_R2.fastq.gz


### Trim reads with bbduk
mkdir ./data/fastq_trimmed
DIR_RAW="./data/fastq"
DIR_TRIM="./data/fastq_trimmed"
#
FQ_FILES=( $(lh $DIR_RAW/*_R[12].fastq.gz | rev | cut -d'/' -f1 | rev | sed -e 's/_R[12].fastq.gz//' | uniq) )
for FQ in ${FQ_FILES[*]}; do
    bash ./scripts/trim_PE_reads.sh $DIR_RAW/${FQ}_R1.fastq.gz \
                                    $DIR_RAW/${FQ}_R2.fastq.gz \
                                    $DIR_TRIM/${FQ}_R1_trim.fastq.gz \
                                    $DIR_TRIM/${FQ}_R2_trim.fastq.gz
done


### STAR index
# STAR index was generated using the following code, where genome fasta file includes the ERCC 
# spike-in sequences, and the GTF annotations were downloaded from ensemble, version GRCz11.98.
STAR \
  --runThreadN 20 \
  --runMode genomeGenerate \
  --genomeDir $STAR_INDEX_DIR \
  --genomeFastaFiles $GENOME_FASTA \
  --sjdbGTFfile $GTF_FILE --sjdbOverhang 100


### Map read with STAR (STAR_2.5.2b)
mkdir ./data/BAMs
SAMPLE_LIST=( $(ls -lh $DIR_TRIM/*_R1_trim.fastq.gz | rev | cut -d'/' -f1 | rev | sed 's/_R1_trim.fastq.gz//' | tr '\n' ' ') )
DIR_TRIM="./data/fastq_trimmed"
DIR_BAM="./data/BAMs"
THREADS=6
for S in ${SAMPLE_LIST[*]}; do
    READS_PAIR_1=$DIR_TRIM/${S}_R1_trim.fastq.gz
    READS_PAIR_2=$DIR_TRIM/${S}_R2_trim.fastq.gz
    FILE_PREFIX=$DIR_BAM/${S}_
    #
    bash ./scripts/map_STAR.sh $STAR_INDEX_DIR $READS_PAIR_1 $READS_PAIR_2 $FILE_PREFIX $THREADS
done

# Rename files and move log files in a separate folder
BAMs=( ./data/BAMs/*_Aligned.out.bam )
for BAM in ${BAMs[*]}; do
    mv $B ${BAM/_Aligned.out.bam/.bam}
done
mkdir ./data/BAMs/LOGs/
mkdir ./data/BAMs/SJ_out/
mkdir ./data/BAMs/Chimeric_out/
mv ./data/BAMs/*Log.final.out ./data/BAMs/LOGs/
mv ./data/BAMs/*Log.out ./data/BAMs/LOGs/
mv ./data/BAMs/*Log.progress.out ./data/BAMs/LOGs/
mv ./data/BAMs/*SJ.out.tab ./data/BAMs/SJ_out/
mv ./data/BAMs/*Chimeric.out.junction ./data/BAMs/Chimeric_out/
mv ./data/BAMs/*Chimeric.out.sam ./data/BAMs/Chimeric_out/

# Sort and index
BAMs=( ./data/BAMs/*_Aligned.out.bam )
THREADS=1
for BAM in ${BAMs[*]}; do
    SAMPLE=$( basename $BAM | sed 's/.bam//')
    sambama_sort.sh $BAM $THREADS
done

###############